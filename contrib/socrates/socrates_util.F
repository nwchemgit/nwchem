c===============================================================c
c                                                               c
c       Socrates - an NWChem module for teaching                c
c       Copyright © 2009 Jeff Hammond                           c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c===============================================================c
c
c     prints shell batch of 2-e AO integrals
c
      subroutine socrates_print_ao2e(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     1                               buf,tol2e)
      implicit none
#include "stdio.fh"
      integer i, j, k, l
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision tol2e
      double precision buf(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
      do i = ilo, ihi
         do j = jlo, jhi
            do k = klo, khi
               do l = llo, lhi
                  if (abs(buf(l,k,j,i)).gt.tol2e) then
                     write(LuOut,1) i,j,k,l,buf(l,k,j,i)
 1                   format(' ao ',4i5,f20.14)
                  endif
               enddo
            enddo
         enddo
      enddo
      call util_flush(LuOut)
      end
c
c     print simple vector of doubles
c
      subroutine socrates_print_dvec(sz,buf)
      implicit none
#include "stdio.fh"
      integer i, sz
      double precision buf(sz)
      write(LuOut,20)
 20   format(/,'------------------------------'/)
      write(LuOut,21)
 21   format('      i       buf(i)')
      write(LuOut,20)
      do i = 1, sz
         write(LuOut,22) i,buf(i)
 22      format(i8,f20.14)
      enddo
      write(LuOut,20)
      call util_flush(LuOut)
      end
c
c     print simple vector of integers
c
      subroutine socrates_print_ivec(sz,buf)
      implicit none
#include "stdio.fh"
      integer i, sz
      integer buf(sz)
      write(LuOut,20)
 20   format(/,'------------------------------'/)
      write(LuOut,21)
 21   format('      i       buf(i)')
      write(LuOut,20)
      do i = 1, sz
         write(LuOut,22) i,buf(i)
 22      format(i8,f20.14)
      enddo
      write(LuOut,20)
      call util_flush(LuOut)
      end

c $Id$
