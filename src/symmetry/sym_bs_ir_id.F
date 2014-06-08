      subroutine sym_bas_irrep_id(basis, oprint, mbf, njr,
     &                            nbf_ir_id, nbf_ir_sc)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "bas.fh"
#include "geom.fh"
#include "geomP.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer basis             ! [input] basis handle
      logical oprint            ! [input]
      integer mbf               ! [input]
      integer njr               ! [input]
      integer nbf_ir_id(mbf,njr)               ! [output] 
      double precision nbf_ir_sc(mbf,njr)      ! [output] 
c     
c     Return the number of functions per irreducible representation
c     that are present in the specified basis set
c     
      integer geom, nir
      integer icent, jcent, ncent
c     
      integer ang_max
      integer maxf, maxr
      integer l_r, k_r
      integer ishlo, ishhi, ilo, ihi, ish, nbf
      integer class, op, ir, itrace, i
      double precision trace
c     
      integer maxireps, maxops, mcent
      parameter (maxireps = 20, maxops=120)
      double precision c(maxireps)
c     
      integer op_class_index(maxops)
c     
      integer nop,  iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps), molgroup_name
      double precision chars(maxireps*maxireps)
c
      integer itype, nprim, ngen, sphcart, n, irrep_dim, m, indx
      logical do_spherical
c
      double precision sym_trace_bas_op
      external sym_trace_bas_op
 
      integer ncent_unique
      integer centlist(100), jc
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
      if (.not. geom_ncent_unique(geom,ncent_unique))
     $     call errquit('ftx:geom?',0, GEOM_ERR)
      maxf = ((ang_max+1)*(ang_max+2))/2 ! Size of largest cart. shell
      maxr = (ang_max+1)*maxf*maxf ! Size of transformation matrices
      do_spherical = bas_is_spherical(basis)
c     
      if (.not. ma_push_get(MT_DBL, maxr, 'r', l_r, k_r))
     $     call errquit('sym_bas_irreps: ma failed', maxr, MA_ERR)
c
c     Get character table information and class mapping
c     
      call sym_group_name(geom,molgroup_name)
      if (.not. sym_char_table(molgroup_name, nop, nir, iclass_dim,
     $     zir, zclass, chars)) 
     $     call errquit('sym_bas_irreps: no char table',0, GEOM_ERR)
      call sym_op_classify(geom, op_class_index)
c
      call ifill(nir*mbf, 0.0d0, nbf_ir_id, 1)
      call dfill(nir*mbf, 0.0d0, nbf_ir_sc, 1)
c
c     Loop thru operations and get the representation matrices
c     
      do ncent = 1, ncent_unique
         icent = unique_cent(ncent,geom)
         if (.not. bas_ce2cnr(basis, icent, ishlo, ishhi))
     $      call errquit('sym_bas_irreps: bas_ce2cnr',0, BASIS_ERR)
         do ish = ishlo, ishhi
            if (.not. bas_continfo(basis, ish, itype, nprim,
     $         ngen, sphcart)) call errquit
     $         ('sym_bas_irreps: bas_continfo', icent,BASIS_ERR)
            if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $         call errquit('sym_bas_irreps: bas_ce2cnr',0,BASIS_ERR)
            n = (ihi-ilo+1)/ngen ! We're doing the first one then copy
            do m = 1, n
               call dfill(nir, 0.0d0, c, 1)
               c(1) = 1
               mcent = 0
               do op = 1, nop-1
                  call sym_bas_op(geom, op, dbl_mb(k_r), maxf, ang_max,
     &                         do_spherical)
                  class = op_class_index(op)
                  jcent = sym_center_map(geom, icent, op)
                  if (icent.eq.jcent) then
                     indx = itype*maxf*maxf           ! Start of type   
                     indx = indx + maxf*(m-1)+(m-1)   ! Diagonal for m
                     c(class) = c(class) + dbl_mb(k_r+indx)
                  else
                     do jc = 1, mcent
                        if (centlist(jc).eq.jcent) goto 11
                     enddo
                     mcent = mcent + 1
                     centlist(mcent) = jcent
 11                  continue
                  endif
               enddo    ! End of operations applied to function
               do ir = 1, nir
                  irrep_dim = nint(chars((ir-1)*nir+1))
                  trace = ddot(nir,c,1,chars((ir-1)*nir+1),1)/dble(nop)
                  if (trace.gt.0.0d0) then
                     do i=1,ngen
                        nbf_ir_id(ilo+ngen*(i-1)+m-1,ir) = 1
                        nbf_ir_sc(ilo+ngen*(i-1)+m-1,ir) =
     &                            dsqrt(dble(mcent+1))
                     enddo
                  endif
               enddo    ! End of projection onto sym characters
 200           continue
            enddo       ! End of functions in one shell
         enddo          ! End of shell
      enddo             ! End of center

      if (oprint) then
         call util_print_centered(6, 'Symmetry analysis of basis',
     $        18, .true.)
         do i = 1, nbf
           write(6,1) 'Basis function',i
           do ir = 1, nir
            if (nbf_ir_id(i,ir).gt.0) then
               write(6,*) zir(ir),nbf_ir_sc(i,ir)
            endif
1           format(a14,i5,1x,a12,1x,a8)
           enddo
         enddo
         do ir = 1, nir
            write(6,*) 'Irrep',zir(ir)
            do i=1,nbf
               if (nbf_ir_id(i,ir).gt.0) write(6,*) 'Bf',i
            enddo
         enddo
      endif
c     
      if (.not. ma_pop_stack(l_r)) call errquit
     $     ('sym_bas_irreps: memory corrupt?',0, MA_ERR)
c
      end
c $Id$
