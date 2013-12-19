      subroutine dft_printinfo()
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "cdft.fh"

      logical oprint
      integer nbf_per_ir(20)


      oprint = util_print('information', print_low)

      if (ga_nodeid().eq.0 .and. oprint) then
         if (util_print('basis summary', print_default)) then
            if (.not.bas_summary_print(ao_bas_han))
     &           call errquit('dft: basis summary print failed',911,
     &       BASIS_ERR)
         else
            write(6,*)
         endif
         if (oadapt .and. util_print('char table', print_high)) then
            call sym_print_char_table(geom)
         endif
         if (oadapt .and. util_print('basis sym', print_default)) then
            call sym_bas_irreps(ao_bas_han, .true., nbf_per_ir)
         endif
         call util_flush(6)
      endif
c
c     Print out other info as requested by user
c
      if (ga_nodeid() .eq. 0) then
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $        call errquit('dft: geom_print', 0, GEOM_ERR)
            write(6,*)
         endif
         if (util_print('symmetry', print_debug)) then
            call sym_print_all(geom, .true., .true., .true., 
     $           .true., .true.)
            write(6,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(ao_bas_han)) 
     $           call errquit('dft: bas_print', 0, BASIS_ERR)
            write(6,*)
         endif
         if (util_print('basis labels',print_high) .and.
     $        ga_nodeid().eq.0) then
            call bas_print_labels(ao_bas_han)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(ao_bas_han)) 
     $           call errquit('dft:gmp', 0, BASIS_ERR)
            write(6,*)
         endif
         call util_flush(6)
      endif
         return
         end
