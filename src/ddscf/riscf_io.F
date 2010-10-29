#ifdef RISCF_USE_DRA
      subroutine riscf_put_dra ( da_three, scr, indx, nij, Nff, npair, 
     $             ndone, g_three, req )
*
* $Id$
*
      implicit none
#include "errquit.fh"

#include "dra.fh"
#include "msgids.fh"

      integer da_three, indx, nij, Nff, npair, ndone, g_three, req
      double precision scr
      dimension scr(nij)

      integer first, last, gfirst, glast, istat

      if (nij .ne. 0) then
        first = indx
        last = first + nij - 1
        gfirst = indx - ndone
        glast = gfirst + nij - 1
        call ga_put (g_three, gfirst, glast, 1, Nff, scr, nij)
      else
        first = 1000000
        last = 0
      endif
      call ga_igop ( msg_riscf_put_1, first, 1, 'min' )
      call ga_igop ( msg_riscf_put_2, last,  1, 'max' )

      glast = last - ndone
      istat = dra_write_section (.false., g_three, 1, glast, 1, Nff,
     $                               da_three, first,  last, 1, Nff, 
     $                               req )
      if ( istat .ne. 0 ) 
     $     call errquit ('requested end after print', istat ,
     &       UNKNOWN_ERR)
      ndone = last

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine riscf_dra_print ( da )
C     prints contents of da to stdout
C     only for debugging purposes, needs global array of full size

      implicit none

#include "global.fh"
#include "stdio.fh"
#include "dra.fh"

      integer da, istat

      integer ga, req
      logical status
      integer type, dim1, dim2
      character*80 name, file

      istat = dra_inquire(da, type, dim1, dim2, name, file)
*ga:1:0
      status = ga_create(type, dim1, dim2, name, dim1, 0, ga)
      if ( status ) then
        istat = dra_read(ga, da, req)
        call ga_print ( ga )
        status = ga_destroy ( ga )
      else
        call util_print_centered (luout,
     $       'could not allocate GA for writing disk array')
      endif
#endif
*
      end
