       subroutine ccsd_driver( rtdb, geom, basis, movecs_in )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "context.fh"
c
c
c
       integer rtdb                        ! [input] database handle
       integer geom, basis
       character*255 movecs_in
c
c
c
       character*255 title_vecs, basis_vecs
       double precision nuclear_charge, charge
       integer nbf_vecs, nsets
       integer k_occ, l_occ, k_evals, l_evals, g_vecs
       integer k_epair, l_epair
       integer npairs
       integer nbf, nsh, nelec, nact, nocc, nopen, nvir, nocc_correl
       integer occ_frozen, occ_lo, occ_hi
       integer vir_frozen, vir_lo, vir_hi, nmo, nmo2
       integer noper, oseg_lo, oseg_hi, npass, blen, max_noper, nocct
       integer noper_core
       integer g_coul, g_exch, g_int1t, g_epair
       integer i,j,ij
       integer ii
       logical oexch,ocoul,oprint,status,ohalf
       double precision emp2, e2pair
       double precision tpass
       character*7 algo
capr
      integer lnoo,lnov,lnvv,lnoov,lnooo,ndima,lnobj,lnovv
      common/len/lnoo,lnov,lnvv,lnoov,lnooo,ndima,lnobj,lnovv
c
       integer maxit,convi,iprt,irstrt,ncor,ndel,l_mo,k_mo,g_objv,
     &         l_aint,k_aint,l_scr,k_scr
capr
c
c
c
       data g_coul,g_exch,g_int1t/-1,-1,-1/
       data oexch,ocoul/.true.,.true./
       data oprint/.true./
       data ohalf/.true./
capr
       data algo/'twofold'/
capr
       data charge/0.d0/
       double precision zip,one,two
       data zip/0.0d00/,one/1.0d00/,two/2.0d00/
c
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
       status = ma_set_hard_fail(.true.)
       status = ma_set_auto_verify(.true.)
c
c
       status = geom_nuc_charge(geom, nuclear_charge)
       status = status.and.bas_numbf(basis, nbf)
       status = status.and.bas_numcont(basis,nsh)
       if (.not. status)
     $   call errquit('moint_driver: failed to get info', 0)
c
c Reconcile number of electrons and occupied/virtual orbitals
c
       nelec = int(nuclear_charge - charge)
       status = MA_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ)
       call dfill(nbf, zip, dbl_mb(k_occ), 1)
       if (.not. rtdb_get(rtdb, 'scfoccupation', MT_DBL, nbf,
     $      dbl_mb(k_occ))) then
         nocc = nelec/2
         nopen = mod(nelec,2)
         nact = nopen
         do i = 1, nocc
            dbl_mb(i-1+k_occ) = two
         enddo
         do i = 1, nact
           dbl_mb(i-1+nocc+k_occ) = one
         enddo
       else
         nocc = 0
         nact = 0
         do i=1,nbf
           if (dbl_mb(k_occ+i-1).eq.two) then
             nocc = nocc + 1
           else if (dbl_mb(k_occ+i-1).ne.zip) then
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
     $    call errquit('mointdriver: cannot read MO vectors header',0)
       if (nbf.ne.nbf_vecs)
     $    call errquit(
     $    'mointdriver: number of basis fun does not match input',0)
       status = MA_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals)
       if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $                       dbl_mb(k_evals), g_vecs)) call
     $      errquit('mointdriver: cannot read MO vectors',0)
c
c Get frozen core & virtual offsets
c
       if (.not.rtdb_get(rtdb,'fourindex:occ_frozen',MT_INT,1,
     $      occ_frozen)) occ_frozen = 0
       if (.not.rtdb_get(rtdb,'fourindex:vir_frozen',MT_INT,1,
     $      vir_frozen)) vir_frozen = 0
c
c Set index ranges
c
       occ_lo = occ_frozen + 1
       occ_hi = nocc
       nocc_correl = occ_hi - occ_lo + 1
c$$$       vir_lo = nocc + 1
c$$$       vir_hi = nbf - vir_frozen
       vir_lo = 1
       vir_hi = nbf
       nvir = nbf - nocc
       nmo2 = vir_hi - vir_lo + 1
       noper = occ_hi - occ_lo + 1
c
c Get Block-length
c
       if (.not.rtdb_get(rtdb,'fourindex:block',MT_INT,1,blen))
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
       if ((ocoul).and.(oexch)) nocct = 2*nocct
       max_noper = moints_numoper( basis, algo, nbf, nocc, nmo2,
     $                             blen, .false.)
C       MAX_NOPER = 10
capr change gt to ge
       if (max_noper.ge.nocct) then
         noper =  (nocc_correl*(nocc_correl+1))/2
       else
capr
         call errquit('cannot hold J and K in core ',0)
capr
         if ((ocoul).and.(oexch)) then
           noper = max_noper/2
         else
           noper = max_noper
         endif
       endif
       noper_core = noper
       if ((ocoul).and.(oexch)) noper_core = noper*2
capr
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
      lnoov=lnoo*nvir
      lnooo=lnoo*nocc
      lnovv=nocc*lnvv
      ndima=lnoo*nocc*(nocc+1)/2
      lnobj=max(lnoov+lnoov+lnoov,nbf*lnoo)
      lnobj=max(lnoov+lnoov+lnooo,lnobj)
      if (.not.ga_create(MT_DBL,lnobj,nvir+1,'objv',lnobj,1,g_objv))
     &     call errquit('ga_create g_objv failed',0)
capr
c
c Allocate operator matrices
c
       if (ocoul)
     $      g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb oper')
*rak     $      g_coul = ga_create_Jblocked(noper,nbf,'Coulomb oper')
       if (oexch) then
         if ((algo.eq.'twofold').and.(.not.(ohalf))) then
           g_exch = ga_create_JKblocked(noper,nbf,nmo2,'X oper')
*rak           g_exch = ga_create_Kblocked(noper,nbf,nmo2,'X oper')
         else
           g_exch = ga_create_JKblocked(noper,nbf,nbf,'X oper')
*rak           g_exch = ga_create_Kblocked(noper,nbf,nbf,'X oper')
           if (.not.ga_create(MT_DBL,lnoov+lnooo,nbf+1,'int1t',
     &         lnoov+lnooo,1,g_int1t))
     &         call errquit('ga_create g_int1t failed',0)
         endif
       endif
c
c Allocate pair energy array
c
       if (.not.ga_create(MT_DBL,1,((nocc_correl*(nocc_correl+1))/2),
     $      'epair',1,1,g_epair))
     $      call errquit('cannot create epair array',0)
c
c Print info
c
       if (ga_nodeid().eq.0) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Four-Index Transformation',
     $                            30,.true.)
         write(6,966) nbf,nsh,nocc,nocc_correl,blen,occ_lo,occ_hi,
     $                vir_lo,vir_hi,movecs_in, noper_core
 966     format(10x,'Number of basis functions:',10x,i5,/,
     $        10x,'Number of shells:',19x,i5,/,
     $        10x,'Number of occupied orbitals:',8x,i5,/,
     $        10x,'Number of occ. correlated orbitals:',1x,i5,/,
     $        10x,'Block length:',23x,i5,/,
     $        10x,'Superscript MO index range:',6x,i3,1x,'-',i3,/,
     $        10x,'Subscript MO index range:',8x,i3,1x,'-',i3,/,
     $        10x,'MO coefficients read from:',2x,a25,/,
     $        10x,'Number of operator matrices in core:',i5)
         call util_flush(6)
         if (ohalf) then
           write(6,957) 
 957       format(10x,'Half-transformed integrals produced')
         else
           write(6,958)
 958       format(10x,'Half-transformed integrals not produced')
         endif
       endif
c
c Clear 4-index statistics
c
       call moints_stats_clear()

c
c Loop over segments of occupied orbitals
c
       if (ga_nodeid().eq.0) write(6,*)
       npass = 0
       oseg_hi = occ_lo - 1
 22    continue
       tpass = tcgtime()
       npass = npass + 1
       oseg_lo = oseg_hi + 1
       oseg_hi = moints_occrange( noper, occ_lo, oseg_lo, occ_hi )
c
c  J & K build routines
c
       if (algo.eq.'sixfold') then
         call moints_build_6x(basis,  occ_lo, oseg_lo, oseg_hi,
     $                        vir_lo, vir_hi, g_vecs,
     $                        g_coul, ocoul, g_exch, oexch,
     $                        blen )
       else
         call moints_build_2xz( basis,  ohalf, occ_lo, oseg_lo, oseg_hi,
     $                          vir_lo, vir_hi, g_vecs,
     $                          g_coul, ocoul, g_exch, oexch, g_int1t,
     $                          blen )
       endif
c
c
c Evaluate and store pair energies
c
capr insert oexch
       if (oexch) then
        call moints_e2pair( nocc, occ_lo, oseg_lo, oseg_hi, nvir,
     $                      nbf, g_exch, dbl_mb(k_evals), g_epair )
       endif
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
C       call moints_stats_print(algo)
c
c Evaluate mp2 correction from pair energies
c
       if ((oexch).and.ga_nodeid().eq.0) then
         npairs = (nocc_correl*(nocc_correl+1))/2
         status = ma_push_get(MT_DBL, npairs, 'Pair energies',
     $                        l_epair, k_epair)
         ii = 0
         emp2 = zip
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
c Print out pair energies
c 
c$$$         write(6,997)
c$$$ 997     format(//,15x,'Pair Energies',/)
c$$$         write(6,998) (dbl_mb(k_epair+ii),ii=0,npairs-1)
c$$$ 998     format(8f10.5)
c
c Print out total correlation correction
c
         write(6,999) emp2
 999     format(//,15x,31('-'),
     $          /,15x,'MP2 Correction:',f16.10,
     $          /,15x,31('-'),/)
         status = ma_pop_stack(l_epair)
       endif
       call ga_sync()
c
      status = ga_destroy(g_epair)
capr
      iprt=0
      status=status.and.ma_push_get(MT_DBL,ndima,'aint',l_aint,k_aint)
      status=status.and.ma_push_get(MT_DBL,lnobj,'scr',l_scr,k_scr)
      if (.not.status)call errquit('error before ccsd_sortmo',0)
      call ccsd_sortmo(nocc,nvir,nbf,dbl_mb(k_aint),ndima,
     $                 dbl_mb(k_scr),lnobj,g_objv,g_coul,
     $                 g_exch,iprt)
      status=status.and.ga_destroy(g_coul)
      status=status.and.ga_destroy(g_exch)
      status=status.and.ma_pop_stack(l_scr)
      if (.not.status) call errquit('error after ccsd_sortmo',0)
c
      maxit=20
      convi=5
      irstrt=0
      ncor=0
      ndel=0
      call rddcom()
      status=ma_push_get(MT_DBL,(nbf*nbf),'movecs',l_mo,k_mo)
      if (.not.status)call errquit('movecs: failed to alloc', l_mo)
      call ga_get(g_vecs,1,nbf,1,nbf,dbl_mb(k_mo),nbf)
      call aoccsd(basis,ncor,nocc,nvir,ndel,nbf,maxit,convi,iprt,
     &            dbl_mb(k_mo),dbl_mb(k_evals),dbl_mb(k_aint),
     &            blen,g_objv,g_int1t)
      status = ma_pop_stack(l_mo)
      status = status.and.ma_pop_stack(l_aint)
      status=status.and.ga_destroy(g_int1t)
      status=status.and.ga_destroy(g_objv)
      if (.not.status)call errquit('deallocate after aoccsd',0)
capr
capr       if (ocoul) status = ga_destroy(g_coul)
capr       if (oexch) status = ga_destroy(g_exch)
      status = ma_pop_stack(l_evals)
c
c Cleanup
c
      call schwarz_tidy()
      status = status.and.ma_pop_stack(l_occ)
      if (.not.status) call errquit('moints: cannot pop stack?',0)
      status = ga_destroy(g_vecs)
      if (ga_nodeid().eq.0) then
        call MA_summarize_allocated_blocks()
      endif
      call int_terminate()
      call ga_sync()
c
      return
      end
