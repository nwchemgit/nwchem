      subroutine ga_conjg(g_a, g_b)
c   Takes the complex conjugate of g_a and store it in g_b
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer g_a, g_b
c
      integer dtypea, da1, da2
      integer dtypeb, db1, db2
      integer iloa, ihia, jloa, jhia
      integer ilob, ihib, jlob, jhib
      integer iA, lda
      integer iB, ldb
c
      call ga_inquire(g_a, dtypea, da1, da2)
      call ga_inquire(g_b, dtypeb, db1, db2)
      if(dtypea .ne. mt_dcpl .or. dtypeb .ne. mt_dcpl) then
        call ga_error('ga_conj: arrays must be complex')
      end if
      call ga_zero(g_b)
c
      call ga_distribution(g_a, ga_nodeid(), iloa, ihia, jloa, jhia)
      call ga_distribution(g_b, ga_nodeid(), ilob, ihib, jlob, jhib)
c      write(luout,*) "nodeid:", ga_nodeid()
      call ga_access(g_a, iloa, ihia, jloa, jhia, iA, lda)
      call ga_access(g_b, ilob, ihib, jlob, jhib, iB, ldb)
c      write(luout,*) "access"
c      write(luout,*) "lo:", ilo, jlo
c      write(luout,*) "hi:", ihi, jhi
      call conjGA(dcpl_mb(iA), ihia-iloa+1, jhia-jloa+1,
     $            dcpl_mb(iB), ihib-ilob+1, jhib-jlob+1)
c      write(luout,*) "conjugated"
      call ga_release_update(g_a, iloa, ihia, jloa, jhia)
      call ga_release_update(g_b, ilob, ihib, jlob, jhib)
c      write(luout,*) "released"
      call ga_sync()

      end subroutine ga_conjg

      subroutine conjGA(A, arows, acols, B, brows, bcols)
      implicit none
#include "stdio.fh"
c
      integer arows, acols
      complex*16 A(arows, acols)
      integer brows, bcols
      complex*16 B(brows, bcols)
c      write(luout,*) A
      B = conjg(A)
      end subroutine conjGA
