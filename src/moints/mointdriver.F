       subroutine moint_driver( rtdb, geom, basis, movecs_in )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
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
       logical status
       double precision nuclear_charge, charge
       integer k_occ, l_occ, k_evals, l_evals, g_vecs
       integer l_epair, k_epair
       integer i, nbf, nsh, nelec, nact, nocc, nopen, nbf_vecs
       integer mo_indx_lo, mo_indx_hi, nsets, nmo
       integer mo_indx2_lo, mo_indx2_hi, nmo2
       integer noper, nvir
       integer g_coul, g_exch
       logical oexch,ocoul,oprint
       double precision emp2,chunk_factor
       data g_coul,g_exch/-1,-1/
       data oexch,ocoul/.true.,.true./
       data oprint/.true./
c
c
c
       logical movecs_read_header, movecs_read, movecs_write
       external movecs_read_header, movecs_read, movecs_write
       integer ga_create_atom_blocked,ga_create_moblocked
       external ga_create_atom_blocked,ga_create_moblocked
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
c
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
c
c Find MO vector filename in the database.
c
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
c
c
       mo_indx_lo = 1
       mo_indx_hi = nocc
       mo_indx2_lo = nocc + 1
       mo_indx2_hi = nbf
C       mo_indx2_lo = 1
C       mo_indx2_hi = nbf
       nvir = nbf - nocc
       nmo2 = mo_indx2_hi - mo_indx2_lo + 1
       noper = mo_indx_hi - mo_indx_lo + 1
c
c
c
       call int_init(1, basis)
       call int_normalize(basis)
       call schwarz_init(geom,basis)
c
c
c
c
c Get tasksize factor
c
       if (.not.rtdb_get(rtdb,'chunkfactor',MT_DBL,1,chunk_factor))
     $      chunk_factor = 1.d0
c
c
c
       if (ga_nodeid().eq.0) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Four-Index Transformation',
     $                            30,.true.)
         write(6,966) nbf,nsh,nocc,mo_indx_lo,mo_indx_hi,
     $                mo_indx2_lo,mo_indx2_hi,movecs_in,chunk_factor
 966     format(10x,'Number of basis functions:',10x,i5,/,
     $        10x,'Number of shells:',19x,i5,/,
     $        10x,'Number of occupied orbitals:',8x,i5,/,
     $        10x,'Superscript MO index range:',6x,i3,1x,'-',i3,/,
     $        10x,'Subscript MO index range:',8x,i3,1x,'-',i3,/,
     $        10x,'MO coefficients read from:',2x,a25,/,
     $        10x,'Chunksize factor:',14x,f10.2)
         call util_flush(6)
       endif
c
c
c
       if (ocoul)
     $      g_coul = ga_create_moblocked(noper,nmo2,'Coulomb operator')
       if (oexch)
     $      g_exch = ga_create_moblocked(noper,nmo2,'exchange operator')
c
c
c       
       call moints_build_JK(rtdb,basis,
     $                      mo_indx_lo,mo_indx_hi,
     $                      mo_indx2_lo,mo_indx2_hi,
     $                      g_vecs,
     $                      g_coul,ocoul,
     $                      g_exch,oexch,
     $                      oprint,
     $                      chunk_factor)

C       call moints_print_opermatrix(noper,nmo2,g_coul)
c
c
c
       if (oexch) then
         status = ma_push_get(MT_DBL,(noper*(noper+1))/2,'pair e2',
     $                                l_epair,k_epair)
         call moints_mp2_eval(nocc,nvir,dbl_mb(k_evals),g_exch,
     $                      emp2,dbl_mb(k_epair))
         status = ma_pop_stack(l_epair)
       endif
c
c
c
       call schwarz_tidy()
       if (ocoul) status = ga_destroy(g_coul)
       if (oexch) status = ga_destroy(g_exch)
       status = ma_pop_stack(l_evals)
       status = status.and.ma_pop_stack(l_occ)
       if (.not.status) call errquit('moints: cannot pop stack?',0)
c
c
       return
       end


