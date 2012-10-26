       subroutine GA_NormF(g_a, value)
C$Id$
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer g_a
       double precision value
c
C     Frobenius norm of a global array
c
       integer ilo, ihi, jlo, jhi, idim, jdim, ld
       integer ind
       Double Precision Work
       Double Precision dlange
       External dlange
c
c     Each process first determines max of data local to it
c
       call ga_sync()
       value = 0.d0
       call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
       if (ihi.gt.0 .and. ihi.ge.ilo) then
          idim = ihi - ilo + 1
          jdim = jhi - jlo + 1
          call ga_access(g_a, ilo, ihi, jlo, jhi, ind, ld)
          Value = dlange('F', IDim, JDim, Dbl_MB(Ind), LD, Work)
          call ga_release(g_a, ilo, ihi, jlo, jhi)
       endif
c
c     Combine partial results
C     NOTE:  dlange is careful about overflows, but what's below is not!
c
       Value = Value * Value
       call ga_sync()
       call ga_dgop(Msg_NormF, value, 1, '+')
       call ga_sync()
       Value = Sqrt( Value )
c
       Return
       end
c
c-----------------------------------------------------------------------
c
       subroutine NGA_NormF(g_a, value)
C$Id$
c
c      This subroutine essentially does the same as ga_normf but it
c      calculates the Frobenius norm of a global array with any number
c      of dimensions.
c
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer g_a
       double precision value
       integer jlo(ga_max_dim), jhi(ga_max_dim), jdim(ga_max_dim)
       integer ld(ga_max_dim), ndim, nelm, i
       integer inc(ga_max_dim)
       integer ind,ind1,ind2,ind3,ind4,ind5,ind6,ind7
       integer i7,i6,i5,i4,i3,i2,i1
       Double Precision Work
       Double Precision dlange
       External dlange
c
c     Each process first determines max of data local to it
c
       if (ga_max_dim.ne.7)
     +   call errquit("nga_normf: ga_max_dim is not 7",0,UERR)
c
       call ga_sync()
       value = 0.d0
       ndim = ga_ndim(g_a)
       call nga_distribution(g_a, ga_nodeid(), jlo, jhi)
       do i = ndim+1, ga_max_dim
         jlo(i) = 1
         jhi(i) = 1
       enddo
       nelm = 1
       do i = 1, ga_max_dim
         jdim(i) = max(0, jhi(i) - jlo(i) + 1)
         nelm = nelm * jdim(i)
       enddo
       value = 0.0d0
       if (nelm.gt.0) then
          call nga_access(g_a, jlo, jhi, ind, ld)
          do i = ndim, ga_max_dim
            ld(i) = 1
          enddo
          do i = 2, ga_max_dim
            ld(i) = ld(i)*ld(i-1)
          enddo
          inc(1) = 1
          do i = 2, ga_max_dim
            inc(i) = ld(i-1)
          enddo
          do i7 = 1, jdim(7)
            ind7 = ind+inc(7)*(i7-1)
            do i6 = 1, jdim(6)
              ind6 = ind7+inc(6)*(i6-1)
              do i5 = 1, jdim(5)
                ind5 = ind6+inc(5)*(i5-1)
                do i4 = 1, jdim(4)
                  ind4 = ind5+inc(4)*(i4-1)
                  do i3 = 1, jdim(3)
                    ind3 = ind4+inc(3)*(i3-1)
                    do i2 = 1, jdim(2)
                      ind2 = ind3+inc(2)*(i2-1)
                      do i1 = 1, jdim(1)
                        ind1 = ind2+inc(1)*(i1-1)
                        value = value + dbl_mb(ind1)*dbl_mb(ind1)
                      enddo ! i1
                    enddo ! i2
                  enddo ! i3
                enddo ! i4
              enddo ! i5
            enddo ! i6
          enddo ! i7
          call nga_release(g_a, jlo, jhi)
       endif
c
c     Combine partial results
c
       call ga_sync()
       call ga_dgop(Msg_NormF, value, 1, '+')
       call ga_sync()
       Value = Sqrt( Value )
c
       Return
       end
c
c-----------------------------------------------------------------------
c
       subroutine NGA_NormF_patch(g_a, ilo, ihi, value)
C$Id$
c
c      This subroutine essentially does the same as ga_normf but it
c      calculates the Frobenius norm of a global array patch with any
c      number of dimensions.
c
c      This is a collective operation.
c
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer g_a
       integer ilo(ga_max_dim), ihi(ga_max_dim)
       double precision value
       integer jlo(ga_max_dim), jhi(ga_max_dim), jdim(ga_max_dim)
       integer ld(ga_max_dim), ndim, nelm, i
       integer inc(ga_max_dim)
       integer ind,ind1,ind2,ind3,ind4,ind5,ind6,ind7
       integer i7,i6,i5,i4,i3,i2,i1
c
c      Each process first determines max of data local to it
c
       if (ga_max_dim.ne.7)
     +   call errquit("nga_normf: ga_max_dim is not 7",0,UERR)
c
       call ga_sync()
       value = 0.d0
       ndim = ga_ndim(g_a)
       call nga_distribution(g_a, ga_nodeid(), jlo, jhi)
       do i = ndim+1, ga_max_dim
         jlo(i) = 1
         jhi(i) = 1
       enddo
c
c      Next establish which part of the local block is within the
c      patch specified on input.
c
       do i = 1, ndim
         jlo(i) = max(jlo(i),ilo(i))
         jhi(i) = min(jhi(i),ihi(i))
       enddo
c
c      Calculate the number of elements of patch that are local to this
c      processor.
c
       nelm = 1
       do i = 1, ga_max_dim
         jdim(i) = max(0, jhi(i) - jlo(i) + 1)
         nelm = nelm * jdim(i)
       enddo
       if (nelm.gt.0) then
          call nga_access(g_a, jlo, jhi, ind, ld)
          do i = ndim, ga_max_dim
            ld(i) = 1
          enddo
          do i = 2, ga_max_dim
            ld(i) = ld(i)*ld(i-1)
          enddo
          inc(1) = 1
          do i = 2, ga_max_dim
            inc(i) = ld(i-1)
          enddo
          do i7 = 1, jdim(7)
            ind7 = ind+inc(7)*(i7-1)
            do i6 = 1, jdim(6)
              ind6 = ind7+inc(6)*(i6-1)
              do i5 = 1, jdim(5)
                ind5 = ind6+inc(5)*(i5-1)
                do i4 = 1, jdim(4)
                  ind4 = ind5+inc(4)*(i4-1)
                  do i3 = 1, jdim(3)
                    ind3 = ind4+inc(3)*(i3-1)
                    do i2 = 1, jdim(2)
                      ind2 = ind3+inc(2)*(i2-1)
                      do i1 = 1, jdim(1)
                        ind1 = ind2+inc(1)*(i1-1)
                        value = value + dbl_mb(ind1)*dbl_mb(ind1)
                      enddo ! i1
                    enddo ! i2
                  enddo ! i3
                enddo ! i4
              enddo ! i5
            enddo ! i6
          enddo ! i7
          call nga_release(g_a, jlo, jhi)
       endif
c
c      Combine partial results
c
       call ga_sync()
       call ga_dgop(Msg_NormF, value, 1, '+')
       call ga_sync()
       Value = Sqrt( Value )
c
       Return
       end
c
c-----------------------------------------------------------------------
c
       double precision function dNormF_patch(d_a, ndim, ilo, ihi, ld)
C$Id$
c
c      This function calculates the Frobenius norm of an array patch
c      with up to 7 dimensions.
c
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       double precision d_a(*)
       integer ndim
       integer ilo(ndim), ihi(ndim), ld(ndim)
       double precision value
       integer jlo(ga_max_dim), jhi(ga_max_dim), jdim(ga_max_dim)
       integer jld(ga_max_dim), nelm, i
       integer inc(ga_max_dim)
       integer ind,ind1,ind2,ind3,ind4,ind5,ind6,ind7
       integer i7,i6,i5,i4,i3,i2,i1
c
c      Each process first determines max of data local to it
c
       if (ga_max_dim.ne.7)
     +   call errquit("dnormf: ga_max_dim is not 7",0,UERR)
       if (ndim.gt.ga_max_dim)
     +   call errquit("dnormf: ndim too large",ndim,UERR)
c
       value = 0.d0
       do i = ndim+1, ga_max_dim
         jlo(i) = 1
         jhi(i) = 1
         jld(i) = 1
       enddo
c
c      Complete the data on which patch we are interested in.
c
       do i = 1, ndim
         jlo(i) = ilo(i)
         jhi(i) = ihi(i)
         jld(i) = ld(i)
       enddo
c
c      Calculate the number of elements of patch that are local to this
c      processor.
c
       nelm = 1
       do i = 1, ga_max_dim
         jdim(i) = max(0, jhi(i) - jlo(i) + 1)
         nelm = nelm * jdim(i)
       enddo
       if (nelm.gt.0) then
          do i = 2, ga_max_dim
            jld(i) = jld(i)*jld(i-1)
          enddo
          inc(1) = 1
          do i = 2, ga_max_dim
            inc(i) = jld(i-1)
          enddo
          ind = 1
          do i7 = 1, jdim(7)
            ind7 = ind+inc(7)*(i7-1)
            do i6 = 1, jdim(6)
              ind6 = ind7+inc(6)*(i6-1)
              do i5 = 1, jdim(5)
                ind5 = ind6+inc(5)*(i5-1)
                do i4 = 1, jdim(4)
                  ind4 = ind5+inc(4)*(i4-1)
                  do i3 = 1, jdim(3)
                    ind3 = ind4+inc(3)*(i3-1)
                    do i2 = 1, jdim(2)
                      ind2 = ind3+inc(2)*(i2-1)
                      do i1 = 1, jdim(1)
                        ind1 = ind2+inc(1)*(i1-1)
                        value = value + d_a(ind1)*d_a(ind1)
                      enddo ! i1
                    enddo ! i2
                  enddo ! i3
                enddo ! i4
              enddo ! i5
            enddo ! i6
          enddo ! i7
       endif
c
c      Compute the final answer
c
       dnormf_patch = Sqrt( Value )
c
       Return
       end
