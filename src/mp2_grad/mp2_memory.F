      subroutine mp2_memory(basis,geom,rtdb,tol2e,oskel,
     $     nbf,noa,nva,nob,nvb,
     $     batch_size,max_ma,max_ga,max_file,task)
*
* $Id: mp2_memory.F,v 1.14 1997-06-28 05:15:50 d3g681 Exp $
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer nbf,basis,rtdb,geom
      integer noa,nva
      double precision tol2e
      logical oskel
      integer nob,nvb
      integer max_ma
      integer max_ga
      double precision max_file
      integer nnbf, nnbf_sparse
      integer nva_local, nsh, ilo, ihi, jlo, jhi, ish,jsh
      character*(*)task
      double precision k_file_size,t_file_size,p_file_size
      integer shdim(1),shlo(1),shpairs(1),shpairslocal(1)
      integer map(1024)  ! No. of processors
      integer lengd,lenscrd
      integer nshpairlocal,nshpair,batch_size,shmax
      integer nbfpairlocal,nbfpair
c
c  Hacked to impose 2 GB limit on a single file.
c     
      nnbf=nbf*(nbf+1)/2
      nva_local=(nva-1)/ga_nnodes()+1
c     
c     Count no. of non-zero basis function pairs in same way as 
c     moints-semi to estimate the file size
c
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('mp2_memory: bas_numcont failed', 0)
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      nnbf_sparse = 0
      do ish = 1, nsh
         if (.not. bas_cn2bfr(basis,ish,ilo,ihi)) call errquit
     $        ('mp2_memory: bas_cntobfr failed', ish)
         do jsh = 1, ish
            if (.not. bas_cn2bfr(basis,jsh,jlo,jhi)) call errquit
     $           ('mp2_memory: bas_cntobfr failed', jsh)
            if (schwarz_shell(ish,jsh)*schwarz_max().ge.tol2e) then
               if (ish.ne. jsh) then
                  nnbf_sparse = nnbf_sparse+(ihi-ilo+1)*(jhi-jlo+1)
               else
                  nnbf_sparse = nnbf_sparse+(ihi-ilo+1)*(ihi-ilo+2)/2
               endif
            endif
         enddo
      enddo
      if (nnbf_sparse .gt. nnbf) call errquit('mp2mem: confused',0)
      call schwarz_tidy()
      call int_terminate()
c
      k_file_size=dble(nnbf_sparse)*batch_size*nva_local
      if (k_file_size .gt. 2.4d8) k_file_size = 1d100 
      if(task.eq.'energy')then
         max_file=k_file_size
      else
         call int_init(rtdb,1,basis)
         call schwarz_init(geom,basis)
         call mp2_backt_info(basis,tol2e,oskel,
     $        .false.,nshpair,nshpairlocal,
     $        nbfpair,nbfpairlocal,shmax,
     $        shdim,shlo,shpairs,shpairslocal,map)
         call int_terminate()
         call schwarz_tidy()
         t_file_size=dble((nob*nvb+noa*nva))*batch_size*nva_local
         p_file_size=dble(nbfpairlocal)*nbf*batch_size
         if (t_file_size .gt. 2.4d8) t_file_size = 1d100 
         if (p_file_size .gt. 2.4d8) p_file_size = 1d100 
         max_file=max(k_file_size+t_file_size,t_file_size+p_file_size)
      endif
c     
      if(task.eq.'energy')then
         max_ga=2*nbf*nbf+noa*nva+nob*nvb
      else
         max_ga=6*nbf*nbf+max(noa*nva+nob*nvb,nbfpair*nva)
      endif
      max_ga = max_ga / ga_nnodes()
c     
      max_ma=2*nbf*nbf+((nva+nvb)*8+nnbf+1)
     $     /ma_sizeof(mt_dbl,1,mt_int)
      if(task.eq.'energy')then
         max_ma=max_ma+2*nbf*nbf+max(noa*nva,nob*nvb)
      else
         call intd_init(rtdb,1,basis)
         call intb_mem_2e4c(lengd, lenscrd) ! blocking algorithm
         lengd = 15**4          ! MUST MATCH MP_NONSEP
         call int_terminate
c
         max_ma=max_ma+max(
     $        (2*nbf*nbf+2*(noa*nva+nob*nvb)),
     $        (4*nbf*nbf+3*(noa*nva+nob*nvb)),
     $        (nshpair*2+nshpairlocal*2+lengd*4) / 
     $            ma_sizeof(mt_dbl,1,mt_int)
     $        +nbf*nbf+2*max(nbf*nbf,nbf*shmax*shmax)
     $        +shmax*shmax*batch_size*nbf+shmax*shmax*nbf
     $        +batch_size*nbf+lenscrd+shmax**4+12*lengd)
      endif
c     
      max_ga = max_ga * 1.1d0    ! 10% for safety
      max_ma = max_ma * 1.1d0
      max_file = max_file * 1.1d0
c     
      end
