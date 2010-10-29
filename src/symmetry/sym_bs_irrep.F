      subroutine sym_bas_irreps(basis, oprint, nbf_per_ir)
C$Id$
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer basis             ! [input] basis handle
      logical oprint            ! [input]
      integer nbf_per_ir(*)     ! [output] 
c     
c     Return the number of functions per irreducible representation
c     that are present in the specified basis set
c     
      integer geom
      integer icent, jcent, ncent
      integer nir
c     
      integer ang_max
      integer maxf, maxr
      integer l_r, k_r
      integer ishlo, ishhi, ilo, ihi, ish, nbf
      integer class, op, ir, itrace, i
      double precision trace
c     
      integer maxireps, maxops
      parameter (maxireps = 20, maxops=120)
      double precision c(maxireps)
c     
      integer op_class_index(maxops)
c     
      integer nop,  iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps), group_name
      double precision chars(maxireps*maxireps)
c
      integer itype, nprim, ngen, sphcart, n, irrep_dim
      logical do_spherical
c
      double precision sym_trace_bas_op
      external sym_trace_bas_op
c
c     Get basis and geom info
c     
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_bas_irreps: no geometry from basis?', basis, BASIS_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('sym_bas_irreps: bas_numbf?', basis, BASIS_ERR)
      if (.not. bas_high_angular(basis, ang_max))
     $     call errquit('sym_bas_irreps: bas_high_angular', basis,
     &       BASIS_ERR)
      maxf = ((ang_max+1)*(ang_max+2))/2 ! Size of largest cart. shell
      maxr = (ang_max+1)*maxf*maxf ! Size of transformation matrices
      do_spherical = bas_is_spherical(basis)
c     
      if (.not. ma_push_get(MT_DBL, maxr, 'r', l_r, k_r))
     $     call errquit('sym_bas_irreps: ma failed', maxr, MA_ERR)
c
      if (.not. geom_ncent(geom, ncent)) call errquit
     $     ('sym_bas_irreps: geom_ncent ?', geom, GEOM_ERR)
c     
c     Get character table information and class mapping
c     
      call sym_group_name(geom,group_name)
      if (.not. sym_char_table(group_name, nop, nir, iclass_dim,
     $     zir, zclass, chars)) 
     $     call errquit('sym_bas_irreps: no char table',0, GEOM_ERR)
      call sym_op_classify(geom, op_class_index)
c
      call dfill(nir, 0.0d0, c, 1)
      c(1) = nbf                ! Have to do indentity manually
c
c     Loop thru operations and get the representation matrices
c     
      do op = 1, nop-1          ! Don't include identity
         call sym_bas_op
     &      (geom, op, dbl_mb(k_r), maxf, ang_max, do_spherical)
         class = op_class_index(op)
c
*         write(6,*) ' OP = ', op, class
c     
         do icent = 1, ncent
c
*            write(6,*) '      CENT = ', icent
c
            if (.not. bas_ce2cnr(basis, icent, ishlo, ishhi))
     $           call errquit('sym_bas_irreps: bas_ce2cnr',0, BASIS_ERR)
            jcent = sym_center_map(geom, icent, op)
            if (icent .eq. jcent) then
c     
c     Center maps into iself ... compute the trace of the matrix
c     
               do ish = ishlo, ishhi
                  if (.not. bas_continfo(basis, ish, itype, nprim,
     $                 ngen, sphcart)) call errquit
     $                 ('sym_bas_irreps: bas_continfo', icent,
     &       BASIS_ERR)
                  if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $                 call errquit('sym_bas_irreps: bas_ce2cnr',0,
     &       BASIS_ERR)
                  n = (ihi - ilo + 1)
                  trace = 0.0d0
                  if (itype .eq. -1) then
                     trace = trace + sym_trace_bas_op(dbl_mb(k_r),
     $                    ang_max, maxf, 1, 0)
                     trace = trace + sym_trace_bas_op(dbl_mb(k_r),
     $                    ang_max, maxf, 3, 1)
                     ngen = 1
                  else
                     n = n / ngen
                     trace = trace + sym_trace_bas_op(dbl_mb(k_r),
     $                    ang_max, maxf, n, itype)
                  endif
c
*                  write(6,*) '         ish ', ish, itype, trace, ngen
c     
                  c(class) = c(class) + trace*ngen
               enddo
            endif
         enddo
      enddo
*c
*      write(6,*) ' NBF = ', nbf
*c     
c
c     Now project onto the characters of the irreps.  All would
c     be simple except for groups where the characters are
c     naturally complex and are combined here as separate irreps
c     whose characters are the real and imaginary parts of one
c     of the complex conjugate pair of irreps.  This leads to the
c     second of this pair having a ZERO character for the identity.
c     Also, if this character is used in a projection operator then
c     the operator is NOT idempotent.  Ugh.   However, twice the
c     character of the first 'irrep' will serve to provide a good 
c     projection operator provided it is faked as if of unit dimension.
c     
      call ifill(nir, 0, nbf_per_ir, 1)
      n = 0
      do ir = 1, nir
         irrep_dim = nint(chars((ir-1)*nir+1))
         if (irrep_dim .eq. 0) then
c
c     This is the wierd case noted above.  Simply reuse the trace
c     from the previous irrep ... it should be the same since the
c     basis set posseses the point group symmetry
c
            continue
         else
            trace = ddot(nir, c, 1, chars((ir-1)*nir+1), 1)/dble(nop)
            itrace = nint(trace)
            if (abs(trace-dble(itrace)).gt.1d-6.or.trace.lt.-1d-6) then
               write(6,*) ' sym_bas_irreps ', (c(i),i=1,nir)
               write(6,*) ' sym_bas_irreps ', trace
               write(6,*) ' sym_bas_irreps ', ir
               call errquit('sym_bas_irreps: bad character', 0,
     &       INPUT_ERR)
            endif
            itrace = itrace*irrep_dim
         endif
         nbf_per_ir(ir) = itrace
         n = n + nbf_per_ir(ir)
      enddo
c
c     Check that we have all of the functions
c
      if (n .ne. nbf) then
         write(6,*) ' sym_bas_irreps ', (c(i),i=1,nir)
         write(6,*) ' sym_bas_irreps ', n, nbf
         call errquit('sym_bas_irreps: lost some functions', 0,
     &       BASIS_ERR)
      endif
c     
      if (oprint) then
         call util_print_centered(6, 'Symmetry analysis of basis',
     $        18, .true.)
         write(6,*)
         do ir = 1, nir
            write(6,1) zir(ir), nbf_per_ir(ir)
         enddo
 1       format(8x, a8, i5)
         write(6,*)
      endif
c
      if (.not. ma_pop_stack(l_r)) call errquit
     $     ('sym_bas_irreps: memory corrupt?',0, MA_ERR)
c
      end
