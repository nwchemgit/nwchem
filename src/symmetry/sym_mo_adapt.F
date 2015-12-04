      subroutine sym_movecs_adapt(basis, thresh, g_vecs, irs, nmixed)
C$Id$
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "global.fh"
#include "util.fh"
#include "msgids.fh"
      integer basis             ! [input]
      double precision thresh   ! [input]
      integer g_vecs            ! [input]
      integer irs(*)            ! [output]
      integer nmixed            ! [output]
c     
c     Symmetry adapt the molecular orbitals in ga_vecs, returning
c     in irs(i) the number of the irreducible representation of
c     the ith molecular orbital.  In nmixed is returned the 
c     number of input molecular orbitals that were symmetry
c     contaminated greater than thresh.  An MO is contaminated if
c     it contains two or more irreps. with coefficients greater than
c     thresh.
c     
c     !! The input MOs are assumed linearly independent.
c     
c     !! It may be necessary to call ga_orthog to reorthogonalize
c     !! the output vectors.
c     
c     !! In the presence of large contamination 
c     !! the output MOs may be exactly linearly dependent.
c     
      integer op 
      integer nop               ! No. of ops including indentity
      integer nir
      integer geom, type, nbf, nvec
      integer g_t, g_q, g_s, l_v, k_v, l_q, k_q, l_t, k_t, me, nproc
      integer l_c, k_c
      integer nvec_per_proc, irrep_dim
      integer class, i, j, ioff, ir, ii, i1, i2
      double precision v1, v2, trace, scale
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      double precision small_neg
      parameter (small_neg=-1.0d-6) ! For negative projection sanity test
c     
      integer maxireps, maxops
      parameter (maxireps=20, maxops=120)
      integer op_class_index(maxops)
      integer iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      character*8 group_name
      double precision v(maxireps), vnorm
      integer idamax
      external idamax
c
      logical odebug            ! True if debugging
      logical owarn             ! True if to print warning messages
c
      odebug = util_print('sym_movecs_adapt_debug',print_never)
      owarn = util_print('sym_movecs_adapt_warn',print_debug)
c
      if (odebug) then
         write(6,*) ' sym_movecs_adapt input orbitals'
         call ga_print(g_vecs)
      endif
c     
c     Get char table and group class info
c     
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_movecs_adapt: no geometry from basis?', basis,
     &       BASIS_ERR)
      call sym_op_classify(geom, op_class_index)
      call sym_group_name(geom,group_name)
      if (.not. sym_char_table(group_name, nop, nir, iclass_dim,
     $     zir, zclass, chars)) 
     $     call errquit('sym_print_char_table: no char table',0,
     &       INPUT_ERR)
c
      if (odebug) then
         do i = 0, ga_nodeid()
            if (i .eq. ga_nodeid()) then
               write(6,*) ' NODE ', i
               call sym_print_all
     $              (geom, .true., .true., .true., .true., .true.)
               call util_flush(6)
            endif
            call ga_sync()
         enddo
      endif
c     
c     Allocate global and local memory
c
      call ga_inquire(g_vecs, type, nbf, nvec)
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nvec, 'sym:tmp',nbf,0,g_t))
     $     call errquit('sym_movecs_adapt: unable to alloc g_t',
     $     nbf*nvec, GA_ERR)
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nvec, 'sym:tmp',nbf,0,g_q))
     $     call errquit('sym_movecs_adapt: unable to alloc g_q',
     $     nbf*nvec, GA_ERR)
c
      g_s = ga_create_atom_blocked(geom, basis,'sym: S')

*ga:1:0
*      if (.not. ga_create(mt_dbl, nbf, nbf, 'sym:S',0,0,g_s))
*     $     call errquit('sym_movecs_adapt: unable to alloc g_s',0)
      if (.not. ma_push_get(mt_dbl, nbf, 'sym:v', l_v, k_v))
     $     call errquit('sym_movecs_adapt: unable to alloc v',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf, 'sym:q', l_q, k_q))
     $     call errquit('sym_movecs_adapt: unable to alloc q',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf, 'sym:t', l_t, k_t))
     $     call errquit('sym_movecs_adapt: unable to alloc t',0, MA_ERR)
      me = ga_nodeid()
      nproc = ga_nnodes()
      nvec_per_proc = (nvec-1)/nproc + 1
      if (.not. ma_push_get(mt_dbl, nir*nvec_per_proc, 'sym:c', 
     $     l_c, k_c))
     $     call errquit('sym_movecs_adapt: unable to alloc c',0, MA_ERR)
c     
      call ga_zero(g_s)
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)
c     
      nmixed = 0
      call ifill(nvec, 0, irs, 1) ! Zeroed for later igop()
      call dfill(nvec_per_proc*nir, 0.0d0, dbl_mb(k_c), 1)
c     
      do class = 1, nir
         if (class.eq.1) then
            call ga_copy(g_vecs, g_t)
         else
            call ga_sync()
            do i = me+1, nvec, nproc
               call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_v), nbf)
               call dfill(nbf, 0.0d0, dbl_mb(k_t), 1)
               do op = 1, nop-1   
                  if (class .eq. op_class_index(op)) then
                     call sym_movecs_apply_op(basis, op, dbl_mb(k_v), 
     $                    dbl_mb(k_t))
                  endif
               enddo
               call ga_put(g_t, 1, nbf, i, i, dbl_mb(k_t), nbf)
            enddo
         endif
c     
c     g_q = Overlap * sum(OP in class) OP*g_vecs
c     
         call ga_dgemm('n','n',nbf,nvec,nbf,1.0d0,g_s, g_t, 0.0d0, g_q)
c     
c     Add into character for given mo and class
c     
         do i = me+1, nvec, nproc
            ii = (i-1)/nproc
            call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_v), nbf)
            call ga_get(g_q,    1, nbf, i, i, dbl_mb(k_q), nbf)
            ioff = k_c + ii*nir + class - 1
            trace = ddot(nbf, dbl_mb(k_v), 1, dbl_mb(k_q), 1)
            dbl_mb(ioff) = dbl_mb(ioff) + trace
         enddo
      enddo
c
      if (odebug) then
         write(6,*) ' Weighted characters of orbitals '
         call output(dbl_mb(k_c), 1, nir, 1, nvec, nir, nvec, 1)
      endif
c     
c     Now dbl_mb(k_c) contains the full character for all local MOs
c     weighted by the number of elements in each class.
c     For each MO project onto the characters of the irreps and
c     determine the dominant irrep and if there is mixing.
c
c     This would be simple if it were not for complex conjugate
c     pairs of irreps that are stored here with one having the
c     real pieces of the characters and the other the imaginary.
c     This leads to the second having a zero character for the
c     identify.  We fudge our way out of this by projecting
c     with twice the real component.
c     
      call ga_sync()
      do i = me+1, nvec, nproc
         ii = (i-1)/nproc
         ioff = k_c + ii*nir
c     
c     Project onto characters
c     
         do ir = 1, nir
            irrep_dim = nint(chars((ir-1)*nir+1))
            if (irrep_dim .eq. 0) then
c
c     This is the wierd case noted above.  Simply pretend that
c     this irrep does not exist and double the previous entry
c     since it will be off by a factor of two ?
c
               v(ir) = 0.0d0
               v(ir-1) = v(ir-1)*2.0d0
            else
               v(ir) = ddot(nir,dbl_mb(ioff),1,chars((ir-1)*nir+1),1)
               scale = dble(irrep_dim) / dble(nop)
               v(ir) = v(ir) * scale
            endif
         enddo
c
         if (odebug) then
            vnorm = 0.0d0
            do ir = 1, nir
               vnorm = vnorm + v(ir)
            enddo
            write(6,*) i, (v(ir),ir=1,nir)
            write(6,*) vnorm
         endif
c
c     Find max component and sanity checks
c     
         do ir = 1, nir
            if (v(ir) .lt. small_neg) then
               write(6,1) i, (v(j),j=1,nir)
 1             format(' sym_movecs_adapt: orbital ',i5,
     $              ' negative proj. '/1x,1p,20d10.2)
               call util_flush(6)
               call errquit('sym_movecs_adapt: negative proj',0,
     &       UNKNOWN_ERR)
            endif
         enddo
c
         i1 = idamax(nir, v, 1)
         v1 = v(i1)
         if (v1 .lt. thresh) call errquit
     $        ('sym_movecs_adapt: no significant component?',i1,
     &       UNKNOWN_ERR)
c     
c     Look for mixing
c     
         v2 = 0.0d0
         if (nir .gt. 1) then
            v(i1) = 0.0d0
            i2 = idamax(nir, v, 1)
            v2 = v(i2)
            v(i1) = v1
         endif
         if (v2 .gt. thresh) then
            nmixed = nmixed + 1
            if (owarn) then
               write(6,11) i, (v(ir),ir=1,nir)
 11            format(' sym_movecs_adapt: orbital ',i5,
     $              ' was contaminated'/1x,1p,20d10.2)
               call util_flush(6)
            endif
         endif
         if (abs(v2-v1) .lt. 0.01d0) then
c
c     Two (or more) irreps are strongly occupied. On the assumption
c     that only two are involved and they are consecutive just
c     assign them to symmetries based on odd/even ... if the assumption
c     is wrong then linear dependence will result.
c
*     call errquit('sym_movecs_adapt: unable to assign MO'//
*     $          ' irrep - try restarting from another geometry',0)
c          if (mod(i,2).eq.0) i1 = i2
       endif
c     
c     Have found the irrep
c
         irs(i) = i1
         ir = i1
c     
c     Only need to project if it was mixed ... but check to high
c     precision (1d-10) not against the user tolerance
c
         if (v2 .gt. 1d-20) then ! Note that V2 is the square of the coeff
c
c     Apply the projection operator ... again, note that
c     in the case of complex conjugate representations we
c     need a factor of two snuck in.
c     
            call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_v), nbf)
c     
c     Identity
c     
           call dcopy(nbf, dbl_mb(k_v), 1, dbl_mb(k_q), 1)
           call dscal(nbf, chars((ir-1)*nir+1), dbl_mb(k_q), 1)
c     
           do op = 1, nop-1
              class = op_class_index(op)
              call dfill(nbf, 0.0d0, dbl_mb(k_t), 1)
              call sym_movecs_apply_op(basis, op, dbl_mb(k_v),
     $             dbl_mb(k_t))
              call daxpy(nbf, chars((ir-1)*nir+class), dbl_mb(k_t), 1,
     $             dbl_mb(k_q), 1)
           enddo
           irrep_dim = nint(chars((ir-1)*nir+1))
           scale = dble(irrep_dim)/dble(nop)
c
c     Handle the weird case of real components of complex chars
c     (Note that we never project with the pure imaginary character)
c
           if (ir .ne. nir) then
              if (chars(ir*nir + 1).eq.0.0d0) scale = scale*2.0d0
           endif
           call dscal(nbf, scale, dbl_mb(k_q), 1)
c       
           if (odebug) then
              write(6,2) i, ir
 2            format('     Input         Output    mo=',i5,', ir=',i2)
              do j = 0,nbf-1
                 write(6,3) dbl_mb(k_v+j), dbl_mb(k_q+j)
 3               format(1x,f13.6,1x,f13.6)
              enddo
           endif
c     
           call ga_put(g_vecs, 1, nbf, i, i, dbl_mb(k_q), nbf)
         endif
c     
      enddo
c     
      call ga_sync
      call ga_igop(msg_sym_1, nmixed, 1, '+')
      call ga_igop(msg_sym_2, irs, nvec, '+')
      call ga_sync
c     
      if (.not. ga_destroy(g_s)) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_q)) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_t)) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
      if (.not. ma_pop_stack(l_c)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_t)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_q)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_v)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
c     
      end
