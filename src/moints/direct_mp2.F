       subroutine direct_mp2(rtdb)
C$Id: direct_mp2.F,v 1.4 1995-02-02 23:18:38 d3g681 Exp $
       implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c
c
c
       integer rtdb                        ! [input] database handle
c
c
c
       integer geom, basis
       character*255 title
       character*255 movecs_in
       character*255 title_vecs, basis_vecs
       double precision nuclear_charge, charge
       integer nbf_vecs, nsets
       integer k_occ, l_occ, k_evals, l_evals, g_vecs
       integer k_epair, l_epair
       integer npairs
       integer nbf, nsh, nelec, nact, nocc, nopen, nvir, nocc_correl
       integer occ_frozen, occ_lo, occ_hi
       integer vir_frozen, vir_lo, vir_hi, nmo, nmo2
       integer noper, oseg_lo, oseg_hi, npass, blen, nocct
       integer g_coul, g_exch, g_epair
       integer i, j, ij, ii
       logical oprint,status
       double precision emp2, e2pair
       double precision tpass, ttotal
       character*7 algo
c
c
c
       data g_coul,g_exch/-1,-1/
       data oprint/.true./
       data algo/'twofold'/
       data charge/0.d0/
c
c
c
       logical movecs_read_header, movecs_read, movecs_write
       external movecs_read_header, movecs_read, movecs_write
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
       integer ga_create_Jblocked,ga_create_Kblocked
       external ga_create_Jblocked,ga_create_Kblocked
       integer moints_numoper
       external moints_numoper
       integer moints_occrange
       external moints_occrange
c
c
c
       if (.not. context_push('direct_mp2'))
     $      call errquit('direct_mp2: context_push failed',0)
       ttotal = tcgtime()
       status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c Memory settings
c
       status = ma_set_hard_fail(.true.)
       status = ma_set_auto_verify(.true.)
c
c Extract high level info from the data-base setting defaults
c     
       if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $      title = ' '
       if (.not. geom_create(geom, 'geometry'))
     $      call errquit('direct_mp2: geom_create?', 0)
       if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $      call errquit('direct_mp2: no geometry ', 0)
       if (.not. bas_create(basis, 'mo basis'))
     $      call errquit('direct_mp2: bas_create?', 0)
       if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $      call errquit('direct_mp2: no mo basis ', 0)
       if (.not. context_rtdb_match(rtdb, 'input mo vector file',
     $      movecs_in)) movecs_in = 'movecs'
       if (.not. geom_nuc_charge(geom, nuclear_charge))
     $      call errquit('direct_mp2: no nuclear charge?', 0)
       if (.not. bas_numbf(basis, nbf))
     $      call errquit('direct_mp2: cannot get basis info',0)
       if (.not. bas_numcont(basis,nsh))
     $      call errquit('direct_mp2: cannot get basis info',0)
c
c Reconcile number of electrons and occupied/virtual orbitals
c
       nelec = int(nuclear_charge - charge)
       status = MA_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ)
       call dfill(nbf, 0.0d0, dbl_mb(k_occ), 1)
       if (.not. rtdb_get(rtdb, 'scfoccupation', MT_DBL, nbf,
     $      dbl_mb(k_occ))) then
         nocc = nelec/2
         nopen = mod(nelec,2)
         nact = nopen
         do i = 1, nocc
            dbl_mb(i-1+k_occ) = 2.0d0
         enddo
         do i = 1, nact
           dbl_mb(i-1+nocc+k_occ) = 1.d0
         enddo
       else
         nocc = 0
         nact = 0
         do i=1,nbf
           if (dbl_mb(k_occ+i-1).eq.2.d0) then
             nocc = nocc + 1
           else if (dbl_mb(k_occ+i-1).ne.0.d0) then
             nact = nact + 1
           endif
         enddo
       endif
C       if (nact.le.0) return
c
c Find MO vector filename in the database.
c
       if (.not.rtdb_cget(rtdb,'moints:input vectors',1,movecs_in))
     $      then
         if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_in))
     $        movecs_in = 'movecs'
       endif
c
c Read in MO coefficients
c
       g_vecs = ga_create_atom_blocked(geom,basis,'Molecular orbital')
       if (.not. movecs_read_header(movecs_in, title_vecs,
     $    basis_vecs, nbf_vecs, nsets, nmo, 2))
     $    call errquit('direct_mp2: cannot read MO vectors header',0)
       if (nbf.ne.nbf_vecs)
     $    call errquit(
     $    'direct_mp2: number of basis fun does not match input',0)
       status = MA_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals)
       if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $                       dbl_mb(k_evals), g_vecs)) call
     $      errquit('direct_mp2: cannot read MO vectors',0)
c
c Get frozen core & virtual offsets
c
       if (.not.rtdb_get(rtdb,'occ_frozen',MT_INT,1,occ_frozen))
     $      occ_frozen = 0
       if (.not.rtdb_get(rtdb,'vir_frozen',MT_INT,1,vir_frozen))
     $      vir_frozen = 0
c
c Set index ranges
c
       occ_lo = occ_frozen + 1
       occ_hi = nocc
       nocc_correl = occ_hi - occ_lo + 1
       vir_lo = nocc + 1
       vir_hi = nbf - vir_frozen
       nvir = nbf - nocc
       nmo2 = vir_hi - vir_lo + 1
c
c Get Block-length
c
       if (.not.rtdb_get(rtdb,'block',MT_INT,1,blen))
     $      blen = 10
c
c Initialise integrals and screening
c
       call int_init(rtdb,1, basis)
       call int_normalize(basis)
       call schwarz_init(geom,basis)
c
c Maximum number of operator matrices we can allocate
c
       nocct = (nocc_correl*(nocc_correl+1))/2
       noper = moints_numoper( basis, algo, nbf, nocc, nmo2, blen)
       if (noper.gt.nocct) noper = nocct
c
c Allocate operator matrix
c
       if (algo.eq.'twofold') then
         g_exch = ga_create_Kblocked(noper,nbf,nmo2,'X oper')
       else
         g_exch = ga_create_Kblocked(noper,nbf,nbf,'X oper')
       endif
c
c Allocate pair energy array
c
       if (.not.ga_create(MT_DBL,1,((nocc_correl*(nocc_correl+1))/2),
     $      'epair',1,1,g_epair))
     $      call errquit('direct_mp2:cannot create epair array',0)
c
c Print info
c
       if (ga_nodeid().eq.0) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Direct MP2', 30, .true. )
         write(6,966) nbf, nsh, blen, occ_lo, occ_hi,
     $                vir_lo, vir_hi, movecs_in, noper
 966     format(
     $        10x,'Basis functions:',20x,i5,/,
     $        10x,'Shells:',29x,i5,/,
     $        10x,'Block length:',23x,i5,/,
     $        10x,'Active occupied range:',11x,i3,1x,'-',i3,/,
     $        10x,'Active virtual range:',12x,i3,1x,'-',i3,/,
     $        10x,'MO coefficients read from:',2x,a25,/,
     $        10x,'Operator matrices in core:',10x,i5,/)
         call util_flush(6)
       endif
c
c Clear 4-index statistics
c
       call moints_stats_clear()
c
c Loop over segments of occupied orbitals
c
       npass = 0
       oseg_hi = occ_lo - 1
 22    continue
       tpass = tcgtime()
       npass = npass + 1
       oseg_lo = oseg_hi + 1
       oseg_hi = moints_occrange( noper, occ_lo, oseg_lo, occ_hi )
c
c J & K build routines: two algorithm choices
c
       if (algo.eq.'sixfold') then
         call moints_build_6x(basis,  occ_lo, oseg_lo, oseg_hi,
     $                        vir_lo, vir_hi, g_vecs,
     $                        g_coul, .false., g_exch, .true.,
     $                        blen )
       else
         call moints_build_2x(basis,  .false., occ_lo,
     $                        oseg_lo, oseg_hi,
     $                        vir_lo, vir_hi, g_vecs,
     $                        g_coul, .false., g_exch, .true.,
     $                        blen )
       endif
c
c Evaluate and store pair energies
c
       call moints_e2pair( nocc, occ_lo, oseg_lo, oseg_hi, nvir,
     $                     g_exch, dbl_mb(k_evals), g_epair )
c
c
c
       tpass = tcgtime() - tpass
       if (ga_nodeid().eq.0) write(6,772) npass, oseg_lo, oseg_hi, tpass
 772   format(5x,'Pass:',i5,5x,'Index range:',i5,2x,'-',i5,
     $        5x,'Time:',f10.2)
       if (oseg_hi.lt.occ_hi) goto 22
c
c End multipass loop
c
       call moints_stats_print(algo)
c
c Evaluate mp2 correction from pair energies
c
       if (ga_nodeid().eq.0) then
         npairs = (nocc_correl*(nocc_correl+1))/2
         status = ma_push_get(MT_DBL, npairs, 'Pair energies',
     $                        l_epair, k_epair)
         ii = 0
         emp2 = 0.0d0
         do i=occ_lo,nocc
           do j=occ_lo,i
             ij = ((i-occ_frozen)*(i-occ_frozen-1))/2 + j - occ_frozen
             call ga_get(g_epair,1,1,ij,ij,e2pair,1)
             dbl_mb(k_epair+ii) = e2pair
             emp2 = emp2 + e2pair
             ii = ii + 1
           enddo
         enddo
c$$$         write(6,997)
c$$$ 997     format(//,15x,'MP2 Pair Energies',/)
c$$$         write(6,998) (dbl_mb(k_epair+ii),ii=0,npairs-1)
c$$$ 998     format(8f10.5)
         write(6,999) emp2
 999     format(//,15x,31('-'),
     $          /,15x,'MP2 Correction:',f16.10,
     $          /,15x,31('-'),/)
         status = ma_pop_stack(l_epair)
       endif
       call ga_sync()
       ttotal = tcgtime() - ttotal
       if (ga_nodeid().eq.0) write(6,823) ttotal
 823   format('Total MP2 time:',5x,f10.2)
c
c Cleanup
c
       status = ga_destroy(g_epair)
       status = ga_destroy(g_exch)
       status = ma_pop_stack(l_evals)
       call schwarz_tidy()
       status = status.and.ma_pop_stack(l_occ)
       if (.not.status) call errquit('direct_mp2: cannot pop stack?',0)
       status = ga_destroy(g_vecs)
       status = geom_destroy(geom)
       if (.not. context_pop('direct_mp2'))
     $     call errquit('context_pop failed',0)
c
c
c
       return
       end


