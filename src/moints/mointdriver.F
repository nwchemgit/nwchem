       subroutine moint_driver( rtdb, geom, basis, movecs_in )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
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
       integer i, nbf, nelec, nact, nocc, nopen, nbf_vecs
       integer mo_indx_lo, mo_indx_hi, nsets, nmo
       integer mo_indx2_lo, mo_indx2_hi, nmo2
       integer noper
       integer g_coul, g_exch
       double precision emp2
c
c
c
       logical movecs_read_header, movecs_read, movecs_write
       external movecs_read_header, movecs_read, movecs_write
       logical ga_create_atom_blocked
       external ga_create_atom_blocked
       integer ga_create_moblocked
       external ga_create_moblocked
c
c
c
       status = ma_set_hard_fail(.true.)
       status = ma_set_auto_verify(.true.)
c
c
       status = geom_nuc_charge(geom, nuclear_charge)
       status = status.and.bas_numbf(basis, nbf)
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
c
       g_vecs = ga_create_atom_blocked(geom,basis,'Molecular orbital')
       if (.not. movecs_read_header(movecs_in, title_vecs,
     $      basis_vecs, nbf_vecs, nsets, nmo, 2)) return
       if (nbf.ne.nbf_vecs) return
       status = MA_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals)
       if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $                       dbl_mb(k_evals), g_vecs)) call
     $      errquit('mointdriver: cannot read MO vectors',0)
c
c
c
       call int_init(1, basis)
       call int_normalize(basis)
c
c
c
       mo_indx_lo = 1
       mo_indx_hi = nocc
C       mo_indx2_lo = nocc + 1
C       mo_indx2_hi = nbf
       mo_indx2_lo = 1
       mo_indx2_hi = nbf
       nmo2 = mo_indx2_hi - mo_indx2_lo + 1
       noper = mo_indx_hi - mo_indx_lo + 1
c
c
c
       if (ga_nodeid().eq.0) then
         write(6,*)
         call util_print_centered(6,'Four-Index Transformation',
     $                            40,.true.)
         write(6,966) nbf,nocc,mo_indx_lo,mo_indx_hi,
     $                mo_indx2_lo,mo_indx2_hi
 966     format(10x,'Number of basis functions:',10x,i5,/,
     $        10x,'Number of occupied orbitals:',8x,i5,/,
     $        10x,'Superscript MO index range:',6x,i3,1x,'-',i3,/,
     $        10x,'Subscript MO index range:',8x,i3,1x,'-',i3)
         call util_flush(6)
       endif

c
c
c
C       g_coul = ga_create_moblocked(noper,nmo2,'Coulomb operator')
       g_exch = ga_create_moblocked(noper,nmo2,'exchange operator')
       call dbg_mark('allocate operator arrays')
c
c
c       
       call moints_build_JK(rtdb,basis,
     $                      mo_indx_lo,mo_indx_hi,
     $                      mo_indx2_lo,mo_indx2_hi,
     $                      g_vecs,g_coul,.false.,g_exch,.true.)
c
c
c
       status = ma_push_get(MT_DBL,(noper*(noper+1))/2,'pair e2',
     $                      l_epair,k_epair)
       call moints_mp2_eval_a(nbf,noper,dbl_mb(k_evals),g_exch,
     $                        emp2,dbl_mb(k_epair))

C       call moints_mp2_eval(nbf,noper,nmo2,dbl_mb(k_evals),
C     $                      g_exch,emp2,dbl_mb(k_epair))
c
c
c
       status = ga_destroy(g_coul)
       status = ga_destroy(g_exch)
       status = ma_pop_stack(l_epair)
       status = ma_pop_stack(l_evals)
       status = status.and.ma_pop_stack(l_occ)
       if (.not.status) call errquit('moints: cannot pop stack?',0)
c
c
       return
       end


