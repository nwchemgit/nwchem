#ifdef RISCF_USE_DA
      subroutine riscf_put_da ( da_three, scr, indx, nij, Nff, npair, 
     $             ndone, g_three )
C$
      implicit none

#include "da.fh"
#include "msgids.fh"

      integer da_three, indx, nij, Nff, npair, ndone, g_three
      double precision scr
      dimension scr(nij)

      integer first, last, gfirst, glast, req, istat

      if (nij .ne. 0) then
        first = indx
        last = first + nij - 1
        gfirst = indx - ndone
        glast = gfirst + nij - 1
C       call ga_print (g_three)
        call ga_put (g_three, gfirst, glast, 1, Nff, scr, nij)
      else
        first = 1000000
        last = 0
      endif
      call ga_igop ( msg_riscf_put_1, first, 1, 'min' )
      call ga_igop ( msg_riscf_put_2, last,  1, 'max' )

      glast = last - ndone
      istat = da_write_section (.false., g_three, 1, glast, 1, Nff,
     $                                da_three, first,  last, 1, Nff, 
     $                                req )
      if ( istat .ne. 0 ) 
     $     call errquit ('requested end after print', istat )
      ndone = last

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine riscf_da_print ( da )
C     prints contents of da to stdout
C     only for debugging purposes, needs global array of full size

      implicit none

#include "global.fh"
#include "da.fh"

      integer da, istat

      integer ga, req
      logical status
      integer type, dim1, dim2
      character*80 name, file

      istat = da_inquire(da, type, dim1, dim2, name, file)
      status = ga_create (type, dim1, dim2, name, dim1, 1, ga)
      istat = da_read(ga, da, req)
      call ga_print ( ga )

      call errquit ('requested end after print',1)
#endif
      end


