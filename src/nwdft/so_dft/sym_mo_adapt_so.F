      subroutine sym_movecs_adapt_so(basis, thresh, g_vecs, irs, nmixed)
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
      integer g_vecs(2)            ! [input]
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
      integer geom, type, nbf, nvec,nbf1,nvec1
      integer g_t(2), g_q(2), g_s, l_v, k_v, l_q, k_q, l_t, k_t, me, 
     &     nproc
      integer g_s2
      integer l_c, k_c
      integer nvec_per_proc, irrep_dim
      integer class, i, j, ioff, ir, ii, i1, i2, pair
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
      integer l_u,k_u
      logical sym_char_table_so
      external sym_char_table_so
c
      logical odebug            ! True if debugging
      logical owarn             ! True if to print warning messages
c
      odebug = util_print('sym_movecs_adapt_debug',print_never)
      owarn = util_print('sym_movecs_adapt_warn',print_debug)
c
      if (odebug) then
c         write(6,*) ' sym_movecs_adapt_so input orbitals'
c         call ga_print(g_vecs(1))
c         call ga_print(g_vecs(2))
      endif
c     
c     Get char table and group class info
c     
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_movecs_adapt_so: no geometry from basis?', basis,
     &       BASIS_ERR)
      nop = sym_number_ops(geom)
      if (.not. ma_push_get(mt_dcpl, 8*nop, 'sym:u', l_u, k_u))
     $     call errquit('sym_movecs_adapt_so: unable to alloc u',
     $                  0, MA_ERR)
c     
      call sym_op_classify_so(geom, op_class_index, dcpl_mb(k_u))
      call sym_group_name(geom,group_name)
      if (.not. sym_char_table_so(group_name, nop, nir, iclass_dim,
     $     zir, zclass, chars)) 
     $     call errquit('sym_print_char_table: no char table',0,
     &       GEOM_ERR)
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
      call ga_inquire(g_vecs(1), type, nbf, nvec)
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nvec, 'sym:tmp',nbf,0,g_t(1)))
     $     call errquit('sym_movecs_adapt_so: unable to alloc g_t',
     $     nbf*nvec, GA_ERR)
      if (.not. ga_create(mt_dbl, nbf, nvec, 'sym:tmp',nbf,0,g_t(2)))
     $     call errquit('sym_movecs_adapt_so: unable to alloc g_t',
     $     nbf*nvec, GA_ERR)
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nvec, 'sym:tmp',nbf,0,g_q(1)))
     $     call errquit('sym_movecs_adapt_so: unable to alloc g_q',
     $     nbf*nvec, GA_ERR)
      if (.not. ga_create(mt_dbl, nbf, nvec, 'sym:tmp',nbf,0,g_q(2)))
     $     call errquit('sym_movecs_adapt_so: unable to alloc g_q',
     $     nbf*nvec, GA_ERR)

c
      g_s = ga_create_atom_blocked(geom, basis,'sym: S')
      call ga_inquire(g_s, type, nbf1, nvec1)
      if (.not. ga_create(mt_dbl, nbf, nbf, 'sym:S',0,0,g_s2))
     $     call errquit('sym_movecs_adapt_so: unable to alloc g_s2',0,
     &       GA_ERR)
*ga:1:0
*      if (.not. ga_create(mt_dbl, nbf, nbf, 'sym:S',0,0,g_s))
*     $     call errquit('sym_movecs_adapt_so: unable to alloc g_s',0)
      if (.not. ma_push_get(mt_dbl, 2*nbf, 'sym:v', l_v, k_v))
     $     call errquit('sym_movecs_adapt_so: unable to alloc v',
     $                  0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, 2*nbf, 'sym:q', l_q, k_q))
     $     call errquit('sym_movecs_adapt_so: unable to alloc q',
     $                  0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, 2*nbf, 'sym:t', l_t, k_t))
     $     call errquit('sym_movecs_adapt_so: unable to alloc t',
     $                  0, MA_ERR)
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      nvec_per_proc = (nvec-1)/nproc + 1
      if (.not. ma_push_get(mt_dbl, 2*nir*nvec_per_proc, 'sym:c', 
     $     l_c, k_c))
     $     call errquit('sym_movecs_adapt_so: unable to alloc c',
     $                  0, MA_ERR)
c     
      call ga_zero(g_s)
      call ga_zero(g_s2)
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)
      call ga_fock_sf(g_s, g_s2, nbf/2)
c      call ga_print(g_s2)
c      call ga_print(g_s)
c     
      nmixed = 0
      call ifill(nvec, 0, irs, 1) ! Zeroed for later igop()
      call dfill(2*nvec_per_proc*nir, 0.0d0, dbl_mb(k_c), 1)
c     
c      call ga_print(g_vecs(1))
c      call ga_print(g_vecs(2))
      do class = 1, nir
         if (class.eq.1) then
            call ga_copy(g_vecs(1), g_t(1))
            call ga_copy(g_vecs(2), g_t(2))
         else if (zclass(class).eq."R") then
            call ga_copy(g_vecs(1), g_t(1))
            call ga_copy(g_vecs(2), g_t(2))
            call ga_scale(g_t(1), -1.0d0)
            call ga_scale(g_t(2), -1.0d0)
         else
            call ga_sync()
            do i = me+1, nvec, nproc 
               if(class.eq.3.and.i.eq.1)then
c                  write(*,*)"test, g_vecs"
c                  call ga_print(g_vecs(1))
c                  call ga_print(g_vecs(2))
               endif
               call ga_get(g_vecs(1),1,nbf,i,i,dbl_mb(k_v),1)
               call ga_get(g_vecs(2),1,nbf,i,i,dbl_mb(k_v+nbf),1)
               call dfill(2*nbf, 0.0d0, dbl_mb(k_t), 1)  
               do op = 1, 2*nop-2 
                  if (class .eq. op_class_index(op))then
                     if(class.eq.3 .and. i.eq.3)then 
c                        write(*,*)"class=,i=", class, i
                     endif
                     call sym_movecs_apply_op_so(basis, op, nop, nbf,   
     $                    dbl_mb(k_v), dbl_mb(k_t), dcpl_mb(k_u))
                     if(class.eq.3 .and. i.eq.3)then 
c                        write(*,*)"k_v_t, op=", op 
c                        do ii=1,nbf*2
c                           write(*,'(2f10.5)')dbl_mb(k_v+ii-1), 
c     &                          dbl_mb(k_t+ii-1)
c                        enddo
                     endif
                  endif
               enddo
               call ga_put(g_t(1),1,nbf,i,i,dbl_mb(k_t),1)
               call ga_put(g_t(2),1,nbf,i,i,dbl_mb(k_t+nbf),1)
            enddo
         endif
c         call ga_print(g_t(1))
c         call ga_print(g_t(2))
c     
c     g_q = Overlap * sum(OP in class) OP*g_vecs
c     
         call ga_dgemm('n','n',nbf,nvec,nbf,1.0d0,g_s2,g_t(1),
     &        0.0d0, g_q(1))
         call ga_dgemm('n','n',nbf,nvec,nbf,1.0d0,g_s2,g_t(2),
     &        0.0d0, g_q(2))
c         call ga_print(g_q(1))
c         call ga_print(g_q(2))
c     
c     Add into character for given mo and class
c     
         do i = me+1, nvec, nproc
            ii = (i-1)/nproc
            call ga_get(g_vecs(1), 1, nbf, i, i, dbl_mb(k_v), 1)
            call ga_get(g_vecs(2), 1, nbf, i, i, dbl_mb(k_v+nbf), 1)
            call ga_get(g_q(1),    1, nbf, i, i, dbl_mb(k_q), 1)
            call ga_get(g_q(2),    1, nbf, i, i, dbl_mb(k_q+nbf), 1)
            ioff = k_c + ii*nir + class - 1
            trace = ddot(nbf, dbl_mb(k_v), 1, dbl_mb(k_q), 1)
c            if(class .eq. 3 .and. i .eq. 1)
c     &           write(*,*)"i,class,trace", i, class, trace
            trace = trace +  
     &           ddot(nbf, dbl_mb(k_v+nbf), 1, dbl_mb(k_q+nbf), 1)
            dbl_mb(ioff) = dbl_mb(ioff) + trace 
c            if(class .eq. 3 .and. i .eq. 1)
c     &           write(*,*)"i,class,trace", i, class, trace 
         enddo
      enddo
c
      odebug=.true.
      if (odebug) then
         write(6,*) ' Weighted characters of orbitals '
         call output(dbl_mb(k_c), 1, nir, 1, nvec, nir, nvec, 1)
      endif
      write(6,*)"character table"
      call output(chars, 1, nir, 1, nir, nir,nir, 1)
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
               scale = dble(irrep_dim) / dble(2*nop)
               v(ir) = v(ir) * scale
            endif
         enddo
c
         if (odebug) then
            vnorm = 0.0d0
            do ir = 1, nir
               vnorm = vnorm + v(ir)
            enddo
            write(*,*)"mo=", i 
            write(6,*)(v(ir),ir=1,nir)
            write(6,'(f10.5)') vnorm
         endif
c
c     Find max component and sanity checks
c     
         do ir = 1, nir
            if (v(ir) .lt. small_neg) then
               write(6,1) i, (v(j),j=1,nir)
 1             format(' sym_movecs_adapt_so: orbital ',i5,
     $              ' negative proj. '/1x,1p,20d10.2)
               call util_flush(6)
               call errquit('sym_movecs_adapt_so: negative proj',0,
     &       GEOM_ERR)
            endif
         enddo
c
         i1 = idamax(nir, v, 1)
         v1 = v(i1)
         if (v1 .lt. thresh) call errquit
     $        ('sym_movecs_adapt_so: no significant component?', 0,
     &       GEOM_ERR)
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
 11            format(' sym_movecs_adapt_so: orbital ',i5,
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
*     call errquit('sym_movecs_adapt_so: unable to assign MO'//
*     $          ' irrep - try restarting from another geometry',0)
            
            if(mod(i-1,2).eq.0)then 
               pair = i1 
            else if(mod(i-1,2).eq.1)then 
               if(i1 .eq. pair)i1 = i2 
            endif
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
 210        continue 
            call ga_get(g_vecs(1), 1, nbf, i, i, dbl_mb(k_v), 1)
            call ga_get(g_vecs(2), 1, nbf, i, i, dbl_mb(k_v+nbf), 1)
c     
c     Identity
c     
            call dcopy(nbf, dbl_mb(k_v), 1, dbl_mb(k_q), 1)
            call dcopy(nbf, dbl_mb(k_v+nbf), 1, dbl_mb(k_q+nbf), 1)
            call dscal(nbf, chars((ir-1)*nir+1), dbl_mb(k_q), 1)
            call dscal(nbf, chars((ir-1)*nir+1), dbl_mb(k_q+nbf), 1)
c     
c     Rotation by 2pi 
c     
            do class = 1,nir
               if(zclass(class).eq."R")goto 211 
            enddo
 211        continue
            call daxpy(nbf, -chars((ir-1)*nir+class), dbl_mb(k_v), 
     $           1, dbl_mb(k_q), 1)
            call daxpy(nbf, -chars((ir-1)*nir+class), 
     $           dbl_mb(k_v+nbf), 1, dbl_mb(k_q+nbf), 1)
c     
            do op = 1, 2*nop-2
               class = op_class_index(op)
               call ga_get(g_vecs(1), 1, nbf, i, i, dbl_mb(k_v), 1)
               call ga_get(g_vecs(2), 1, nbf, i, i, dbl_mb(k_v+nbf), 1)
               call dfill(2*nbf, 0.0d0, dbl_mb(k_t), 1)
               call sym_movecs_apply_op_so(basis, op, nop, nbf,   
     &              dbl_mb(k_v), dbl_mb(k_t), dcpl_mb(k_u))
               call daxpy(nbf, chars((ir-1)*nir+class), dbl_mb(k_t), 
     $              1, dbl_mb(k_q), 1)
               call daxpy(nbf, chars((ir-1)*nir+class), 
     $              dbl_mb(k_t+nbf), 1, dbl_mb(k_q+nbf), 1)
            enddo
            irrep_dim = nint(chars((ir-1)*nir+1))
            scale = dble(irrep_dim)/dble(2*nop)
c     
c     Handle the weird case of real components of complex chars
c     (Note that we never project with the pure imaginary character)
c     
            if (ir .ne. nir) then
               if (chars(ir*nir + 1).eq.0.0d0) scale = scale*2.0d0
            endif
            call dscal(nbf, scale, dbl_mb(k_q), 1)
            call dscal(nbf, scale, dbl_mb(k_q+nbf), 1)
c     
            if (odebug) then
c               write(6,2) i, ir
 2             format('     Input         Output    mo=',i5,', ir=',i2)
               do j = 0,nbf-1
c                  write(6,3) dbl_mb(k_v+j), dbl_mb(k_q+j)
 3                format(1x,f13.6,1x,f13.6)
               enddo
               do j = 0,nbf-1
c                  write(6,3) dbl_mb(k_v+nbf+j), dbl_mb(k_q+nbf+j)
               enddo
            endif
c            scale = 1.0d0/v(ir)  
c            call dscal(nbf, scale, dbl_mb(k_q), 1)
c            call dscal(nbf, scale, dbl_mb(k_q+nbf), 1)
c            trace = ddot(nbf, dbl_mb(k_v), 1, dbl_mb(k_q), 1)
c            trace = trace + 
c     &           ddot(nbf, dbl_mb(k_v+nbf), 1, dbl_mb(k_q+nbf), 1) 
c            scale=sign(scale,trace)  
c            call dscal(nbf, scale, dbl_mb(k_q), 1)
c            call dscal(nbf, scale, dbl_mb(k_q+nbf), 1) 
c     
            call ga_put(g_vecs(1), 1, nbf, i, i, dbl_mb(k_q), 1)
            call ga_put(g_vecs(2), 1, nbf, i, i, dbl_mb(k_q+nbf), 1)
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
      if (.not. ga_destroy(g_s2)) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_q(1))) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_q(2))) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_t(1))) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_t(2))) call errquit('sym_mos_ad: ga?',0,
     &       GA_ERR)
c
      if (.not. ma_pop_stack(l_c)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_t)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_q)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_v)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_u)) call errquit('sym_mos_ad: ma?',0,
     &       MA_ERR)
c     
      end

      subroutine sym_movecs_apply_op_so(basis, op, nop, nbf2, v, t, u)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
      integer basis, op, nop, opp, nbf, i, j, l_tmp, k_tmp, nbf2, ii
      double precision v(*), t(*), fact
      complex*16 u(*), ufact
c
      nbf = nbf2/2 
      if (.not. ma_push_get(mt_dbl, nbf, 'tmp', l_tmp, k_tmp))
     $     call errquit('sym_movecs_apply_so: unable to alloc tmp',0,
     &       MA_ERR)
      do i = 1, 4 
         opp=mod(op-1,nop-1)+1
         call dfill(nbf, 0.0d0, dbl_mb(k_tmp), 1)
         call sym_movecs_apply_op(basis, opp, 
     &        v((i-1)*nbf+1), dbl_mb(k_tmp))  
         do j = 1, 4 
            ufact=u(4*(op-1)+2*mod(i-1,2)+mod(j-1,2)+1)
            if((i-1)/2 .eq. (j-1)/2)then 
               fact = real(real(ufact))
            else 
               fact = -dimag(ufact)
               if((i-1)/2 .eq. 1)fact = -fact 
            endif
c            fact = 0  
            call daxpy(nbf,fact,dbl_mb(k_tmp),1,t((j-1)*nbf+1),1)
            if(fact.ne.0 .and. (opp .eq. 1 .or. opp .eq. 2))then  
c               write(*,*)"i,j,fact", i,j,fact 
               do ii = 1,44 
c                  write(*,'(3f10.5)')v(ii),dbl_mb(k_tmp+ii-1),t(ii)
               enddo 
            endif
         enddo
      enddo
c
      if (.not. ma_pop_stack(l_tmp)) 
     &     call errquit('sym_movecs_apply: ma?',0,
     &       MA_ERR)
c
      end
      
