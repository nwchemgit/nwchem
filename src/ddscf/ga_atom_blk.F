      integer function ga_create_atom_blocked(geom, basis, name)
C$Id: ga_atom_blk.F,v 1.9 2003-10-21 23:51:25 edo Exp $
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
c
c     Allocate a global array that is distributed so that atom
cs     blocks are not spilt between processors.
c
c     Returns the GA handle.  Hard fail upon an error.  Name is the
c     name passed to ga_create.
      integer geom, basis
      character*(*) name
c
      integer k_map,l_map,nblock(2)
      logical status
      integer nbf,g_a
c
      if(.not.ma_push_get(MT_int,ga_nnodes()*2,
     .        'ga block maps',l_map,k_map))
     .     call errquit('gacratblock: not enough local mem',0,MA_ERR)
      call util_atom_blocked(ga_nnodes(),
     .     geom, basis, int_mb(k_map), nblock)
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ga_create_atom_blocked: bas_nbf ?',0,BASIS_ERR)
#ifdef DEBUG
      if(ga_nodeid().eq.0) then
      write(0,*) nblock(1),' m1 ',int_mb(k_map),int_mb(k_map+1)
      write(0,*) nblock(2),' m2 ',int_mb(k_map+ga_nnodes()),
     ,     int_mb(k_map+ga_nnodes()+1)
      do g_a=1,2*ga_nnodes()
      write(0,*) ' i m ',g_a,int_mb(k_map+g_a-1)
      enddo
      endif
#endif
      status =  ga_create_irreg(MT_DBL, nbf, nbf, name,
     $     int_mb(k_map), nblock(1), 
     $     int_mb(k_map+ga_nnodes()), nblock(2), g_a)
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' ga_create_atom_blocked: ', name
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('ga_create_atom_blocked: ga_create_irreg ',nbf,
     $  GA_ERR)
      endif
c
      ga_create_atom_blocked = g_a
      if (.not.ma_chop_stack(l_map))
     &     call errquit('xc_rep_init: cannot chop stack',0,MA_ERR)
c
      end
      subroutine util_atom_blocked(max_nproc,
     ,     geom, basis, map, nblock)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
      integer geom, basis
      integer max_nproc
      integer map(max_nproc,2)
c
      integer chunk1, chunk2, n_in_proc1, n_in_proc2, nblock(2)
      integer iat, ihi, ilo, min_chunk
      integer natoms, nbf
      parameter (min_chunk = 10)
c
c     min_chunk ... minimum chunksize to avoid excessive fragmentation
c     
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('util_atom_blocked: geom_ncent ?', 0, GEOM_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('util_create_atom_blocked: bas_nbf ?', 0, BASIS_ERR)
c
c     chunk1/2 = mean no. of bf / processor need in each dimension
c
      chunk1 = min_chunk
      chunk2 = min_chunk
 10   continue
c
c     Partition basis functions so all those of a given atom are grouped
c     together and that chunks in each dimension exceed the minimal
c     chunking
c
      n_in_proc1 = chunk1
      n_in_proc2 = chunk2
      nblock(1) = 0
      nblock(2) = 0
      do iat = 1, natoms
         if (.not. bas_ce2bfr(basis, iat, ilo, ihi)) call errquit
     $        ('ga_create_atom_blocked: bas_ce2bfr ?', 0, BASIS_ERR)
         if ((ihi-ilo+1) .gt. 0) then
            if (n_in_proc1 .ge. chunk1) then
               nblock(1) = nblock(1) + 1
               if (nblock(1) .le. max_nproc) map(nblock(1),1) = ilo
               n_in_proc1 = 0
            endif
            n_in_proc1 = n_in_proc1 + (ihi-ilo) + 1
            if (n_in_proc2 .ge. chunk2) then
               nblock(2) = nblock(2) + 1
               if (nblock(2) .le. max_nproc) map(nblock(2),2) = ilo
               n_in_proc2 = 0
            endif
            n_in_proc2 = n_in_proc2 + (ihi-ilo) + 1
         endif
      enddo
      if (nblock(1)*nblock(2) .gt. max_nproc) then
         if (chunk1.eq.nbf .and. chunk2.eq.nbf)
     $        call errquit('ga_create_atom_blocked: chunking', 0,
     &       GA_ERR)
         if (chunk1 .lt. nbf) then
            chunk1 = chunk1 + 1
         endif
         if (chunk2 .lt. nbf) then
            chunk2 = chunk2 + 1
         endif
         goto 10
      endif
#ifdef DEBUG
      if(ga_nodeid().eq.0) then
      write(0,*) ' nblock1 ',nblock(1)
      write(0,*) ' map1 ',(map(iat,1),iat=1,nblock(1))
      write(0,*) ' nblock2 ',nblock(2)
      write(0,*) ' map2 ',(map(iat,2),iat=1,nblock(2))
      endif
#endif
      return
      end
