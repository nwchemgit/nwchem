#define DEBUG_MAP_ 1
      integer function ga_create_atom_blocked(geom, basis, name)
C$Id$
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
      integer ilo,ihi,jlo,jhi
c
      if(.not.ma_push_get(MT_int,ga_nnodes()*2,
     .        'ga block maps',l_map,k_map))
     .     call errquit('gacratblock: not enough local mem',0,MA_ERR)
      call util_atom_blocked(ga_nnodes(),
     .     geom, basis, int_mb(k_map), nblock)
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ga_create_atom_blocked: bas_nbf ?',0,BASIS_ERR)
#ifdef DEBUG_MAP
      if(ga_nodeid().eq.0) then
      write(0,*) nblock(1),' m1 ',(int_mb(k_map+g_a),g_a=0,nblock(1)-1)
      write(0,*) nblock(2),' m2 ',(int_mb(k_map+ga_nnodes()+g_a),
     A     g_a=0,nblock(2)-1)
      endif
#endif
      status =  ga_create_irreg(MT_DBL, nbf, nbf, name,
     $     int_mb(k_map), nblock(1), 
     $     int_mb(k_map+ga_nnodes()), nblock(2), g_a)
#ifdef DEBUG_MAP
      call ga_distribution(g_a,
     .     ga_nodeid(), ilo,ihi,jlo,jhi)
      write(6,123) ga_nodeid(),' g2d distr ilo',ilo,ihi,jlo,jhi
 123  format(i2,a,4i5)
#endif
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
      integer function ga_create3d_atom_blocked(geom, basis, name, dim1)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
c
c     3d version of ga_create_atom_blocked
c     first dim is not touched
c
      integer geom, basis
      integer dim1
      character*(*) name
c
      integer k_map,l_map
      logical status
      integer nbf,g_a
      integer dims(3)
      integer nblock(3)
      integer jj
c
      if(.not.ma_push_get(MT_int,ga_nnodes()*2+1,
     .        'ga block maps',l_map,k_map))
     .     call errquit('gacratblock: not enough local mem',0,MA_ERR)
      call util_atom_blocked(ga_nnodes(),
     .     geom, basis, int_mb(k_map+1), nblock(2))
c     need to compact mapping ie get rid of zeroes
      do g_a=1,nblock(3)
         int_mb(k_map+nblock(2)+g_a)=int_mb(k_map+ga_nnodes()+g_a)
      enddo
#ifdef DEBUG_MAP
      if(ga_nodeid().eq.0) then
      write(0,*) nblock(1),' m1 ',int_mb(k_map+1),int_mb(k_map+2)
      write(0,*) nblock(2),' m2 ',int_mb(k_map+ga_nnodes()+1),
     ,     int_mb(k_map+ga_nnodes()+2)
      do g_a=1,2*ga_nnodes()
      write(0,*) '3d i m ',g_a,int_mb(k_map+g_a)
      enddo
      endif
#endif
      int_mb(k_map)=1
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('ga_create3d_atom_blocked: bas_nbf ?',0,BASIS_ERR)
      dims(1)=dim1
      dims(2)=nbf
      dims(3)=nbf
      nblock(1)=1
#ifdef DEBUG_MAP
      write(6,123) ga_nodeid(),' dim1 ',dim1
      write(6,123) ga_nodeid(),' nblock ',nblock
 123  format(i2,a,4i5)
      do jj=1,nblock(1)*nblock(2)*nblock(3)
         write(6,123) ga_nodeid(),'mapp ',jj,int_mb(k_map+jj-1)
      enddo
#endif
      status =  nga_create_irreg(MT_DBL, 
     N     3, dims, name,
     $     int_mb(k_map), nblock, 
     $     g_a)
      if (.not. status) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' ga_create3d_atom_blocked: ', name
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('ga_create3d_atom_blocked: ga_create_irreg ',nbf,
     $  GA_ERR)
      endif
c
      ga_create3d_atom_blocked = g_a
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
C#ifdef DEBUG
C     if(ga_nodeid().eq.0) then
C     write(0,*) ' nbf ',nbf
C     write(0,*) ' nproc ',max_nproc
C     write(0,*) ' nblock1 ',nblock(1)
C     write(0,*) ' map1 ',(map(iat,1),iat=1,nblock(1))
C     write(0,*) ' nblock2 ',nblock(2)
C     write(0,*) ' map2 ',(map(iat,2),iat=1,nblock(2))
C     endif
C#endif
      return
      end
