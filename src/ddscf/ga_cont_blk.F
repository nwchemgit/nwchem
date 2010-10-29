      integer function ga_create_cont_blocked(basis, name)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
      integer basis
      character*(*) name
c
c     Allocate a global array that is distributed so that 
c     contraction/shell blocks are not spilt between processors.
c
      integer max_nproc
      parameter (max_nproc = 1024)
      integer map1(max_nproc), map2(max_nproc)
      integer chunk1, chunk2, n_in_proc1, n_in_proc2, nblock1, nblock2
      integer icont, ihi, ilo, min_chunk, g_a, nchunk
      integer ncont, nbf
      parameter (min_chunk = 10)
      logical status
c
c     min_chunk ... minimum chunksize to avoid excessive fragmentation
c
      if (.not. bas_numcont(basis, ncont)) call errquit
     $     ('ga_create_cont_blocked: bas_numcont ?', 0, BASIS_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ga_create_cont_blocked: bas_nbf ?', 0, BASIS_ERR)
      if (ga_nnodes().gt.max_nproc) call errquit
     $     ('ga_create_cont_blocked: static dim too small', max_nproc,
     &       GA_ERR)
c
c     chunk1/2 = mean no. of bf / processor need in each dimension
c
      chunk1 = min_chunk
      chunk2 = min_chunk
 10   continue
c
c     Partition basis functions so all those of a given contraction/shell
c     are grouped together and that chunks in each dimension exceed the 
c     minimal chunking
c
      n_in_proc1 = chunk1
      n_in_proc2 = chunk2
      nblock1 = 0
      nblock2 = 0
      do icont = 1, ncont
         if (.not. bas_cn2bfr(basis, icont, ilo, ihi)) call errquit
     $        ('ga_create_cont_blocked: bas_ce2bfr ?', 0, BASIS_ERR)
         if (n_in_proc1 .ge. chunk1) then
            nblock1 = nblock1 + 1
            if (nblock1 .le. max_nproc) map1(nblock1) = ilo
            n_in_proc1 = 0
         endif
         n_in_proc1 = n_in_proc1 + (ihi-ilo) + 1
         if (n_in_proc2 .ge. chunk2) then
            nblock2 = nblock2 + 1
            if (nblock1 .le. max_nproc) map2(nblock2) = ilo
            n_in_proc2 = 0
         endif
         n_in_proc2 = n_in_proc2 + (ihi-ilo) + 1
      enddo
c
      nchunk = nblock1*nblock2
      if (nblock1*nblock2 .gt. ga_nnodes()) then
*         if (ga_nodeid() .eq. 0) then
*            write(6,*) ' P, c1, c2, nc ', ga_nnodes(), chunk1, chunk2, 
*     $           nchunk
*            write(6,*) ' nb1, nb2 ', nblock1, nblock2
*            call ffflush(6)
*         endif
         if (chunk1.eq.nbf .and. chunk2.eq.nbf)
     $        call errquit('ga_create_cont_blocked: chunking', 0,
     &       GA_ERR)
         if (chunk1 .lt. nbf) then
            chunk1 = chunk1 + 1
         endif
         if (chunk2 .lt. nbf) then
            chunk2 = chunk2 + 1
         endif
         goto 10
      endif
c
      status =  ga_create_irreg(MT_DBL, nbf, nbf, name,
     $     map1, nblock1, map2, nblock2, g_a)
      if (.not. status)
     $     call errquit('ga_create_cont_blocked: ga_create_irreg ', 0,
     &       GA_ERR)
c
      ga_create_cont_blocked = g_a
c
      end
