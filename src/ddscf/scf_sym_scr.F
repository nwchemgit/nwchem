      subroutine scf_sym_screen(nclosed, nopen, nvir, irs,
     $                          g_v, ivoff, oprint, maxelt)
      implicit none
#include "global.fh"
      integer nclosed, nopen, nvir ! [input] No. of MOs in each space
      integer irs(*) ! [input] symmetry of each MO
      integer g_v    ! [input/output] rotation vector to screen
      integer ivoff  ! [input] Offset within g_v to start of vector
      logical oprint ! [input] If true print warnings
      double precision maxelt ! [output] Returns maximum forbidden element
c
c     Screen the orbital rotation/gradient vector for symmetry
c     forbidden rotations.  Zero them and return the largest.
c
c     Vector is assumed to hold (nclosed+nopen)*nvir + nclosed*nopen
c
      integer i, j, a, ioff
      double precision test
      double precision buf(10000)  ! temporary while debugging
c
      call ga_sync
c
      maxelt = 0.0d0
      do i = 1+ga_nodeid(),nclosed+nopen,ga_nnodes()
        if ((nclosed+nopen)*nvir .gt. 0) then
          ioff = 1+(i-1)*nvir + ivoff
          call ga_get(g_v, ioff, ioff+nvir-1, 1, 1, buf, 1)
          do a=1,nvir
            test = abs(buf(a))
            if (irs(i).ne.irs(a+nclosed+nopen)) then
               maxelt = max(maxelt,test)
               buf(a) = 0.0d0
               if (oprint .and. test.gt.1d-6) then
                 write(6,1) i, a+nclosed+nopen, irs(i),
     $             irs(a+nclosed+nopen), test
 1               format('scf_sym_screen: symmetry forbidden rotation ',
     $              2i5,4x,2i5,4x,d9.2)
                 call util_flush(6)
               endif
            endif
          enddo
          call ga_put(g_v, ioff, ioff+nvir-1, 1, 1, buf, 1)
        endif
c
        if (i.le.nclosed .and. nopen*nclosed.gt.0) then
          ioff = 1+(nclosed+nopen)*nvir+(i-1)*nopen + ivoff
          call ga_get(g_v, ioff, ioff+nopen-1, 1, 1, buf, 1)
          do a=1,nopen
            test = abs(buf(a))
            if (irs(i).ne.irs(a+nclosed)) then
               maxelt = max(maxelt,test)
               buf(a) = 0.0d0
               if (oprint .and. test.gt.1d-6) then
                 write(6,1) i, a+nclosed, irs(i),
     $             irs(a+nclosed), test
                 call util_flush(6)
               endif
            endif
          enddo
          call ga_put(g_v, ioff, ioff+nopen-1, 1, 1, buf, 1)
        endif
      enddo
c
      call ga_sync
c
      end
