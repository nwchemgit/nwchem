      subroutine mp2_memory(basis,geom,rtdb,tol2e,oskel,
     $     nbf,noa,nva,nob,nvb,
     $     batch_size,max_ma,max_ga,max_file,task)
*
* $Id: mp2_memory.F,v 1.12 1997-04-26 21:57:55 d3g681 Exp $
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf,basis,rtdb,geom
      integer noa,nva
      double precision tol2e
      logical oskel
      integer nob,nvb
      integer max_ma
      integer max_ga
      double precision max_file
      integer nnbf
      integer nva_local
      character*(*)task
      double precision k_file_size,t_file_size,p_file_size
      integer shdim(1),shlo(1),shpairs(1),shpairslocal(1)
      integer map(1024)  ! No. of processors
      integer lengd,idum,lenscr,lenscrd,idum2,leng
      integer nshpairlocal,nshpair,batch_size,shmax
      integer nbfpairlocal,nbfpair
c
c  Hacked to impose 2 GB limit on a single file.
c     
      nnbf=nbf*(nbf+1)/2
      nva_local=(nva-1)/ga_nnodes()+1
c     
      k_file_size=dble(nnbf)*batch_size*nva_local
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
      call int_init(rtdb,1,basis)
      call int_mem_2e4c(leng, lenscr)
      call int_terminate
      lenscr = max(lenscr,leng)
c
      max_ma=2*nbf*nbf+((nva+nvb)*8+nnbf+1)
     $     /ma_sizeof(mt_dbl,1,mt_int)
      if(task.eq.'energy')then
         max_ma=max_ma+2*nbf*nbf+max(noa*nva,nob*nvb)
      else
         call int_init(rtdb,1,basis)
         call intd_init(rtdb,1,basis)
         call int_mem ( idum, lengd, idum2, lenscrd )
         lenscrd = max(lenscrd,lengd)
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
