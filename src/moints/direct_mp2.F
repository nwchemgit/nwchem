       subroutine direct_mp2(rtdb)
C$Id: direct_mp2.F,v 1.13 1995-08-07 16:37:17 rg240 Exp $
       implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "util.fh"
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
       integer nbf, nsh, nelec, nact, nocc, nopen, nocc_correl
       integer occ_frozen, occ_lo, occ_hi
       integer vir_frozen, vir_lo, vir_hi, nmo, nmo2
       integer noper, oseg_lo, oseg_hi, npass, blen, nocct
       integer g_coul, g_exch, g_epair
       integer i, j, ij, ii
       logical oprint, status, omemprint
       logical ohalf
       double precision emp2, e2pair
       double precision tpass, ttotal, scf_energy
       character*7 algo
c
c
c
       data g_coul,g_exch/-1,-1/
       data charge/0.d0/
       data ohalf/.false./
c
c
c
       logical movecs_read_header, movecs_read, movecs_write
       external movecs_read_header, movecs_read, movecs_write
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
       integer ga_create_JKblocked
       external ga_create_JKblocked
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
c Print settings
c
       call util_print_push()
       call util_print_rtdb_load(rtdb,'mp2')
       oprint = util_print('information',print_default)
       omemprint = util_print('memory',print_high)
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
       if (.not. bas_create(basis, 'ao basis'))
     $      call errquit('direct_mp2: bas_create?', 0)
       if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $      call errquit('direct_mp2: no ao basis ', 0)
       if (.not. rtdb_cget(rtdb, 'input mo vector file', 1,
     $      movecs_in)) call util_file_prefix('movecs',movecs_in)
       if (.not. geom_nuc_charge(geom, nuclear_charge))
     $      call errquit('direct_mp2: no nuclear charge?', 0)
       if (.not. bas_numbf(basis, nbf))
     $      call errquit('direct_mp2: cannot get basis info',0)
       if (.not. bas_numcont(basis,nsh))
     $      call errquit('direct_mp2: cannot get basis info',0)
       if (.not. rtdb_get(rtdb, 'scf:energy', MT_DBL, 1, scf_energy))
     $      scf_energy = 0.0d0
c
c Get algorithm choice
c
       if (.not. rtdb_cget(rtdb, 'fourindex:method ', 1, algo))
     $      algo = 'twofold'
c
c Reconcile number of electrons and occupied/virtual orbitals
c
       nelec = int(nuclear_charge - charge)
       status = ma_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ)
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
     $        call util_file_prefix('movecs',movecs_in)
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
       if (.not.rtdb_get(rtdb, 'direct_mp2:occ_frozen', MT_INT, 1,
     $                   occ_frozen)) then
         if (.not.rtdb_get(rtdb, 'fourindex:occ_frozen', MT_INT, 1,
     $                     occ_frozen)) occ_frozen = 0
       endif
       if (.not.rtdb_get(rtdb, 'direct_mp2:vir_frozen', MT_INT, 1,
     $                   occ_frozen)) then
         if (.not.rtdb_get(rtdb, 'fourindex:vir_frozen', MT_INT, 1,
     $                     vir_frozen)) vir_frozen = 0
       endif
c
c Set index ranges
c
       occ_lo = occ_frozen + 1
       occ_hi = nocc
       nocc_correl = occ_hi - occ_lo + 1
       vir_lo = nocc + 1
       vir_hi = nbf - vir_frozen
       nmo2 = vir_hi - vir_lo + 1
c
c Get Block-length
c
       if (.not.rtdb_get(rtdb,'fourindex:block',MT_INT,
     $      1,blen)) blen = 10
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
       noper = moints_numoper( basis, algo, nbf, nocc,
     $                         nmo2, blen, omemprint )
       if (noper.gt.nocct) noper = nocct
c
c Allocate operator matrix
c
       if ((algo.eq.'twofold').and.(.not.(ohalf))) then
         g_exch = ga_create_JKblocked(noper,nbf,nmo2,'X oper')
       else
         g_exch = ga_create_JKblocked(noper,nbf,nbf,'X oper')
       endif
c
c Allocate pair energy array
c
       if (.not.ga_create(MT_DBL,1,((nocc_correl*(nocc_correl+1))/2),
     $      'epair',1,1,g_epair))
     $      call errquit('direct_mp2:cannot create epair array',0)
c
c Dummy loop to count number of passes
c
       npass = 0
       oseg_hi = occ_lo - 1
 322   continue
       npass = npass + 1
       oseg_lo = oseg_hi + 1
       oseg_hi = moints_occrange( noper, occ_lo, oseg_lo, occ_hi )
       if (oseg_hi.lt.occ_hi) goto 322
c
c Print info
c
       if ((ga_nodeid().eq.0).and.oprint) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Direct MP2', 30, .true. )
         write(6,966) nbf, nsh, blen, occ_lo, occ_hi,
     $                vir_lo, vir_hi, movecs_in, noper,npass
 966     format(
     $        10x,'Basis functions:',20x,i5,/,
     $        10x,'Shells:',29x,i5,/,
     $        10x,'Block length:',23x,i5,/,
     $        10x,'Active occupied range:',11x,i3,1x,'-',i3,/,
     $        10x,'Active virtual range:',12x,i3,1x,'-',i3,/,
     $        10x,'MO coefficients read from:',2x,a25,/,
     $        10x,'Operator matrices in core:',10x,i5,/,
     $        10x,'AO passes required:',17x,i5)
         if (algo.eq.'sixfold') then
           write(6,975) 
 975       format(10x,'Sixfold algorithm used')
         else
           write(6,976) 
 976       format(10x,'Twofold algorithm used')
         endif
         write(6,*)
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
         call moints_build_2x(basis, ohalf, occ_lo,
     $                        oseg_lo, oseg_hi,
     $                        vir_lo, vir_hi, g_vecs,
     $                        g_coul, .false., g_exch, .true.,
     $                        blen )
       endif
c
c Evaluate and store pair energies
c
       call moints_e2pair( nocc, occ_lo, oseg_lo, oseg_hi, nmo2,
     $                     nbf, g_exch, dbl_mb(k_evals), g_epair )
c
c Print progress of this pass
c
       tpass = tcgtime() - tpass
       if ((ga_nodeid().eq.0).and.
     $     (util_print('progress reports',print_default)))
     $     write(6,772) npass, oseg_lo, oseg_hi, tpass
 772   format(5x,'Pass:',i5,5x,'Index range:',i5,2x,'-',i5,
     $        5x,'Time:',f10.2)
c
c
c
c End multipass loop
c
       if (oseg_hi.lt.occ_hi) goto 22
c
c Print statistics
c
       if (util_print('statistics',print_high))
     $      call moints_stats_print(algo)
c
c Gather pair energies and evaluate mp2 correction
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
c
c Print total and pair energies
c
         if (util_print('pair energies',print_debug)) then
           write(6,997)
 997       format(//,15x,'MP2 Pair Energies',/)
           write(6,998) (dbl_mb(k_epair+ii),ii=0,npairs-1)
 998       format(8f10.5)
         endif
         if (oprint) write(6,999) scf_energy, emp2, scf_energy+emp2
 999     format(//,10x,43('-'),
     $          /,10x,'SCF energy         ',f24.12,
     $          /,10x,'Correlation energy ',f24.12,
     $          /,10x,'Total MP2 energy   ',f24.12,
     $          /,10x,43('-'),/)
         status = ma_pop_stack(l_epair)
       endif
       call ga_sync()
       ttotal = tcgtime() - ttotal
       if ((ga_nodeid().eq.0).and.(oprint)) write(6,823) ttotal
 823   format(/,10x,'Total MP2 time:',5x,f10.2)
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
c
c Pop context
c
       call util_print_pop()
       if (.not. context_pop('direct_mp2'))
     $     call errquit('context_pop failed',0)
c
c
c
       return
       end


