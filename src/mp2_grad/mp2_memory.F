      subroutine mp2_memory(basis,rtdb,tol2e,oskel,
     $     nbf,noa,nva,nob,nvb,
     $     batch_size,max_ma,max_ga,max_file,task)
*
* $Id: mp2_memory.F,v 1.10 1997-03-27 08:50:07 d3e129 Exp $
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf,basis,rtdb
      integer noa,nva
      double precision tol2e
      logical oskel
      integer nob,nvb
      integer max_ma
      integer max_ga
      integer max_file
      integer nnbf
      integer nva_local
      character*(*)task
      integer k_file_size,t_file_size,p_file_size
      integer shdim,shlo,shpairs,shpairslocal,map
      integer lengd,idum,lenscr,lenscrd,idum2,leng
      integer nshpairlocal,nshpair,batch_size,shmax
      integer nbfpairlocal,nbfpair
c     
      nnbf=nbf*(nbf+1)/2
      nva_local=(nva-1)/ga_nnodes()+1
c     
      k_file_size=nnbf*batch_size*nva_local
      if(task.eq.'energy')then
         max_file=k_file_size
      else
         call mp2_backt_info(basis,tol2e,oskel,
     $        .false.,nshpair,nshpairlocal,
     $        nbfpair,nbfpairlocal,shmax,
     $        shdim,shlo,shpairs,shpairslocal,map)
         t_file_size=(nob*nvb+noa*nva)*batch_size*nva_local
         p_file_size=nbfpairlocal*nbf*batch_size
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
      call int_init(rtdb,1,basis)
      call int_mem_2e4c(leng, lenscr)
      call int_terminate
c
      max_ma=2*nbf*nbf+((nva+nvb)*8+nnbf+1)
     $     /ma_sizeof(mt_dbl,1,mt_int)
      if(task.eq.'energy')then
         max_ma=max_ma+2*nbf*nbf+max(noa*nva,nob*nvb)
      else
         call int_init(rtdb,1,basis)
         call intd_init(rtdb,1,basis)
         call int_mem ( idum, lengd, idum2, lenscrd )
         call int_terminate
c
         max_ma=max_ma+max(
     $        (2*nbf*nbf+2*(noa*nva+nob*nvb)),
     $        (4*nbf*nbf+3*(noa*nva+nob*nvb)),
     $        (nshpair*2+nshpairlocal*2)/ma_sizeof(mt_dbl,1,mt_int)
     $        +nbf*nbf+2*max(nbf*nbf,nbf*shmax*shmax)
     $        +shmax*shmax*batch_size*nbf+shmax*shmax*nbf
     $        +batch_size*nbf+lenscr+lenscrd+leng+lengd)
      endif
c     
      max_ga = max_ga * 1.1d0    ! 10% for safety
      max_ma = max_ma * 1.1d0
      max_file = max_file * 1.1d0
c     
      end
