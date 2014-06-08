      subroutine tddft_grad_irreg_acc(ipol,nmat,nbf,ibf,data,buf,
     &                                ibt,iet,g_g)
      implicit none
c
c     $Id$
c
c     This is an unfortunate routine. The situation is that there is
c     an array 'data' defined as
c
c       data(ipol*nmat,nbf,nbf)
c
c     and its values need to be added onto a global array g_g which is
c     distributed according to chunking (ipol*nmat,-1,-1). But there is
c     scatter operation involved:
c
c       do l = 1, nbf
c         do k = 1, nbf
c           do i = 1, ipol*nmat
c             g_g(i,ibf(k),ibf(l)) = data(i,k,l)
c           enddo
c         enddo
c       enddo
c
c     I know of no GA operation that can deal with this situation. So
c     the best that can be achieved is to analyse ibf for blocks of
c     consequtive indeces, copy those blocks to a buffer and use the
c     normal nga_acc to add them onto the global array.
c
c     This subroutine has been extended to deal with the case that 
c     g_g is actually an array of 2D GA's instead of a 3D GA as well.
c
c     Written by Huub van Dam, November 2006.
c
#include "global.fh"
#include "errquit.fh"
c
c     Input:
c
      integer ipol     ! =1 (restricted), =2 (unrestricted)
      integer nmat     ! the number of matrices
      integer nbf      ! the number of active basis functions
      integer ibf(nbf) ! the table giving the true basis function number
                       ! for every element of compressed list
c
      double precision data(ipol*nmat,nbf,nbf) ! the data to add on
c
c     Input/Output:
c
      integer g_g(*)   ! the global array(s) the data will be added onto
c
c     Workspace:
c
      integer ibt(nbf) ! the beginnings of consequetive chunks
      integer iet(nbf) ! the endings    of consequetive chunks
c
      double precision buf(ipol*nmat*nbf*nbf) ! the buffer for nga_acc
c
c     Local:
c
      integer nc       ! the number of consequetive chunks
      integer i, k, l  ! counters
      integer icl, ick ! counters over consequetive chunks
      integer lenk     ! the number of consequetive functions in the 
                       ! k-index
      integer ioffk    ! the offset for the k-index in the buffer
      integer ioffl    ! the offset for the l-index in the buffer
c
      integer ld(2)    ! the leading dimensions of buf for nga_acc
      integer ilo(3)   ! the lower limits for the block in g_g
      integer ihi(3)   ! the upper limits for the block in g_g
c
      integer ndim     ! the dimension of the global array g_g
c
c     Code:
c
      nc = 1
      ibt(1) = 1
      iet(1) = 1
      do i = 2, nbf
        if (ibf(i).eq.0) then
          goto 10
        else if (ibf(i-1)+1.ne.ibf(i)) then
          nc = nc + 1
          ibt(nc) = i
        endif
        iet(nc) = i
      enddo
 10   continue
c
      ndim = ga_ndim(g_g)
      if (ndim.eq.2) then
c
        do icl = 1, nc
          do ick = 1, nc
c
c           Copy all data for this block to the buffer
c
            do i = 1, ipol*nmat
              do l = ibt(icl), iet(icl)
                ioffl = l-ibt(icl)
                do k = ibt(ick), iet(ick)
                  lenk = iet(ick)-ibt(ick)+1
                  ioffk = k-ibt(ick)+1
                  buf(ioffk+ioffl*lenk) = 
     &                data(i,k,l)
                enddo
              enddo
c
c             Setup the information for the accumulate
c
              ld(1)  = lenk
              ilo(1) = ibf(ibt(ick))
              ihi(1) = ibf(iet(ick))
              ilo(2) = ibf(ibt(icl))
              ihi(2) = ibf(iet(icl))
c
c             Accumulate this block of consequetive data onto the global
c             array
c
              call nga_acc(g_g(i),ilo,ihi,buf,ld,1.0d0)
c
            enddo
c
          enddo
        enddo
c
      else if (ndim.eq.3) then
c
        do icl = 1, nc
          do ick = 1, nc
c
c           Copy all data for this block to the buffer
c
            do l = ibt(icl), iet(icl)
              ioffl = l-ibt(icl)
              do k = ibt(ick), iet(ick)
                lenk = iet(ick)-ibt(ick)+1
                ioffk = k-ibt(ick)
                do i = 1, ipol*nmat
                  buf(i+ioffk*ipol*nmat+ioffl*ipol*nmat*lenk) = 
     &                data(i,k,l)
                enddo
              enddo
            enddo
c
c           Setup the information for the accumulate
c
            ld(1)  = ipol*nmat
            ld(2)  = lenk
            ilo(1) = 1
            ihi(1) = ipol*nmat
            ilo(2) = ibf(ibt(ick))
            ihi(2) = ibf(iet(ick))
            ilo(3) = ibf(ibt(icl))
            ihi(3) = ibf(iet(icl))
c
c           Accumulate this block of consequetive data onto the global
c           array
c
            call nga_acc(g_g,ilo,ihi,buf,ld,1.0d0)
c
          enddo
        enddo
c
      else 
c
        call errquit('tddft_grad_irreg_acc: unexpected dimension',ndim,
     &               UERR)
c
      endif
c
      end
