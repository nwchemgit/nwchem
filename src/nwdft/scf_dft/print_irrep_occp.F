      subroutine print_irrep_occp(basis, nbf_ao, irs, occ, nspin)
*
      implicit none
c
#include "errquit.fh"
#include "bas.fh"
#include "inp.fh"
#include "global.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "cscfps.fh"
      integer basis
      integer nbf_ao
      integer nspin
      integer irs(nbf_ao,nspin)
      double precision occ(nbf_ao,nspin)
c
      integer i, isp
      integer geom
c
      integer maxireps
      parameter (maxireps=20)
      integer nir, nop, iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      double precision totocc(maxireps,2)
      character*8 group_name
      character*255 title
c
      if (ga_nodeid() .eq. 0) then
         if (.not. bas_geom(basis, geom))
     $    call errquit('print_irrep_occp',0, BASIS_ERR)
c
         call sym_group_name(geom, group_name)
         if (.not. sym_char_table(group_name, nop, nir, iclass_dim,
     $        zir, zclass, chars))
     $    call errquit('print_irrep_occp: no char table',0, UNKNOWN_ERR)
c
         if (nir .eq. 1)  return   ! no symmetry
c
         do i = 1, maxireps
           totocc(i,1) = 0.0
           totocc(i,2) = 0.0
         enddo
         if (nspin .eq. 2) then
           do i = 1, nbf_ao
             totocc(irs(i,1),1) = totocc(irs(i,1),1) + occ(i,1)
             totocc(irs(i,2),2) = totocc(irs(i,2),2) + occ(i,2)
           enddo
         else
           do i = 1, nbf_ao
             totocc(irs(i,1),1) = totocc(irs(i,1),1) + 0.5*occ(i,1)
             totocc(irs(i,1),2) = totocc(irs(i,1),2) + 0.5*occ(i,1)
           enddo
         endif
c
         title = 'Occupations of the irreducible representations'
         write (6,*)
         call util_print_centered(6,title, 40, .true.)
         write(6,*)

         write (6,1001) 'irrep   ','   alpha','    beta'
         write (6,1001) '--------','--------','--------'
         do i = 1, nir
           write (6,1000) zir(i), totocc(i,1), totocc(i,2)
         enddo
         write(6,*)
c
      endif
c
 1000 format(21x,a8,2(5x,f8.1))
 1001 format(21x,a8,2(5x,a8))
c
      end
c $Id$
