c
c     Prints the resulting gradients for each of the excited states of interest.
c
      subroutine tddft_grad_print_result(nroots,iroot,nat,ihdl_geom,g_g)
c
c     Format is taken from  nwchem/src/gradients/grad_force.F
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
c
c     Input:
c
      integer nroots        ! the number of roots
      integer iroot(nroots) ! the root numbers 
      integer nat           ! the number of atoms
      integer ihdl_geom     ! the geometry handle
      integer g_g           ! the global array containing the gradients
c
c     Local:
c
      double precision crd(3) ! coordinate of one atom
      integer ilo(3), ihi(3)  ! patch delimiters
      integer ld(3)           ! leading dimensions
      character*16 tag        ! atom label
      double precision q      ! atom charge
      integer j               ! counter
      integer ir              ! counter over roots
      integer iat             ! counter over atoms
      integer l_grd,k_grd     ! handle and index of gradient memory
c
c     Code:
c
      if (ga_nodeid().eq.0) then
        if (.not.ma_push_get(mt_dbl,3*nat,'grad',l_grd,k_grd)) call
     +    errquit('tddft_grad_print_results: could not allocate grad',
     +            0, MA_ERR)
        ld(1)  = 1
        ld(2)  = 3
        ilo(2) = 1
        ilo(3) = 1
        ihi(2) = 3
        ihi(3) = nat
        do ir=1,nroots
          write(luout,1000)iroot(ir),'x','y','z','x','y','z'
 1000     format(/,/,25x,'Root ',i3,/,/,25x,'TDDFT ENERGY GRADIENTS',
     +           /,/,4X,'atom',15X,'coordinates',24X,'gradient',/,6X,
     +           2(1X,(3(10X,A1))))
          ilo(1) = ir
          ihi(1) = ir
          call nga_get(g_g,ilo,ihi,dbl_mb(k_grd),ld)
          do iat=1,nat
            if (.not. geom_cent_get(ihdl_geom,iat,tag,crd,q)) 
     +        call errquit
     +             ('tddft_grad_print_results: geometry corrupt?',0, 
     +              GEOM_ERR)
            write(luout,2000) iat, tag,(crd(j),j=1,3),
     +           (dbl_mb(k_grd+3*(iat-1)+j),j=0,2)
 2000       format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
          enddo
        enddo
        if (.not.ma_pop_stack(l_grd)) call
     +    errquit('tddft_grad_print_results: could not deallocate grad',
     +            0, MA_ERR)
      endif
      end
c $Id$
