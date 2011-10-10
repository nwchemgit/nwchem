c===============================================================c
c                                                               c
c       NCC - a New Coupled-Cluster code for NWChem             c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c               Karol Kowalski                                  c
c               Environmental Molecular Sciences Laboratory     c
c               Pacific Northwest National Laboratory           c
c               karol.kowalski@pnl.gov                          c
c                                                               c
c               Marta Włoch                                     c
c               Department of Chemistry                         c
c               Michigan Technological University               c
c               wloch@mtu.edu                                   c
c                                                               c
c===============================================================c
c
c     prints shell batch of 2-e AO integrals
c
      subroutine ncc_print_ao2e(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     1                           buf,tol2e)
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
                  if (abs(buf(l,k,j,i)).ge.tol2e) then
                     write(6,1) i,j,k,l,buf(l,k,j,i)
 1                   format(' ao ',4i5,f20.14)
                  endif
               enddo
            enddo
         enddo
      enddo
      call util_flush(LuOut)
      end
c
c     print MO vectors
c
      subroutine ncc_print_movecs(nbf,nmo,buf)
      implicit none
#include "stdio.fh"
      integer ao, mo, nbf, nmo
      double precision buf(nbf,nmo)
      write(6,*)
      write(6,11)
 11   format('      mo      ao       movec(ao,mo)')
      do mo = 1, nmo
         do ao = 1, nbf
            write(6,12) mo,ao,buf(ao,mo)
 12         format(2i8,f20.14)
         enddo
      enddo
      write(6,*)
      call util_flush(LuOut)
      end
c
c     print simple vector of doubles
c
      subroutine ncc_print_dvec(nbf,buf)
      implicit none
#include "stdio.fh"
      integer ao, nbf
      double precision buf(nbf)
      write(6,*)
      write(6,21)
 21   format('      ao       vec(ao)')
      do ao = 1, nbf
         write(6,22) ao,buf(ao)
 22      format(i8,f20.14)
      enddo
      write(6,*)
      call util_flush(LuOut)
      end
c
c     print simple vector of integers
c
      subroutine ncc_print_ivec(nbf,buf)
      implicit none
#include "stdio.fh"
      integer ao, nbf
      integer buf(nbf)
      write(6,*)
      write(6,21)
 21   format('      ao       vec(ao)')
      do ao = 1, nbf
         write(6,22) ao,buf(ao)
 22      format(i8,i20)
      enddo
      write(6,*)
      call util_flush(LuOut)
      end
c
c
c
      subroutine put_ao2e_to_buf(ilo,ihi,jlo,jhi,
     1                           klo,khi,llo,lhi,
     2                           nbf,buf_in,buf_out)
      implicit none
      integer i, j, k, l
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi, nbf
      double precision buf_in(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
      double precision buf_out(1:nbf,1:nbf,1:nbf,1:nbf)
      do i = ilo, ihi
         do j = jlo, jhi
            do k = klo, khi
               do l = llo, lhi
                 buf_out(l,k,j,i) = buf_in(l,k,j,i)
               enddo
            enddo
         enddo
      enddo
      end
c $Id$
