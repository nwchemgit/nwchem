      subroutine ga_chol(g_a, g_s,info)
C$Id: ga_chol.F,v 1.4 1995-02-02 23:19:51 d3g681 Exp $
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
      integer ii
      integer g_a               ! input matrix
      integer g_s               ! output matrix
c
c
c     This subroutine performs a submatrix Choleski factorization of a 
c     symmetric positive matrix with a column wrap using a list of processors.
c
c     
c
c     the input matrix is not destroyed
c
c     proc - # of processors
c     pans - panel size in columns
c     n    - problem size
c
c     
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      integer iproc, pan, istart, iend, k, n, npan, proc, myproc
      integer mypanel, myelem, mycol, elems, ipan, ik, info
      logical status
c
      integer hma, adrma        !mapa
c
      integer hs, adrs          !s
c
      integer his, adris
      integer hip, adrip
c
      external  choleski,fmemreq
c
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimS1, dimS2, typeS
      integer ga_chol_msg
c
      integer Elem
      Elem(istart,iend,n)=((iend-istart+1)*(2*n-istart-iend+2))/2
c
c
c***  check environment
c
      myproc = ga_nodeid()
      proc   = ga_nnodes()
c
c***  check GA info
c
      call ga_check_handle(g_a, 'ga_chol: a')
      call ga_check_handle(g_s, 'ga_chol: s')
      call ga_inquire(g_a,  typeA, dimA1, dimA2)
      call ga_inquire(g_s,  typeS, dimS1, dimS2)
c
c
      if (dimA1 .ne. dimA2) then
         call ga_error('ga_chol: can only fact square matrix ', 1)
      else if (dimS1 .ne. dimS2) then
         call ga_error('ga_chol: can only fact square matrix ', 3)
      else if(dimA1.ne.dimS1 )then
         call ga_error('ga_chol: A & S dimension mismatch ', 4)
      endif
c
      n = dimA1
      call ga_sync()
      proc=ga_nnodes()
      if (n .lt. 350) then
         proc = min(8,ga_nnodes())
      else if (n .lt. 700) then
         proc = min(48,ga_nnodes())
      endif
      oactive = myproc .lt. proc
      if (oactive) then
         pan = 1                !panel size
         npan = n/pan
c     
c        allocate map arrays
c
         status = ma_push_get(mt_int,n,'mapa',hma,adrma)
c         
c        determine distribution of a and z
c     
         mypanel = 0            ! number of panels at given processor
         mycol   = 0            ! number of columns at given processor
         myelem  = 0            ! number of elements at given processor
c     
c        allocate first half of the panels - forwards
c
         do ipan = 1, npan/2 
            iproc = mod(ipan -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = iproc
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c        allocate second half of the panels - backwards
c
         ik = 1                 !ik is used to forward number processors
         do ipan =   npan-1, npan/2 +1 , -1
            ik = ik+1
            iproc = mod(ik -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = iproc
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c***     actually, there is one more panel left for the mismatch 
c     
         iproc  = 0             !processor 0 gets the mismatch panel
         istart = (npan - 1)*pan + 1
         iend   = n
         do k = istart, iend
            int_mb(-1+adrma+k) = iproc
         enddo
         mycol = mypanel * pan
         if(iproc .eq. myproc)then
            myelem = myelem + elem(istart,iend,n)
            mycol  = mycol + iend - istart + 1	
         endif
         elems=myelem
c
c***     allocate a and s arrays
c     
         status = ma_push_get(mt_dbl,elems,'s',hs,adrs)
c
c***     load the fock (A) matrix to (S)
c
         
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. myproc) then
               myelem = n-k+1
               call ga_get(g_a,k,n,k,k,dbl_mb(adrs+istart),1)
               istart = istart + myelem
            endif
         enddo
c     
c***     allocate some memory 
c


         if(.not. ma_push_get(mt_int,(6*n),'iwork',his,adris))
     &        call ga_error('ga_chol allocator',1)
         if(.not. ma_push_get(mt_dbl,(6*n+1),'work',hip,adrip))
     &        call ga_error('ga_chol allocator',3)
c     
      endif
c
c
      call ga_sync()            ! To reduce message buffering problems
c
c***  calling PEIGS parallel choleski
c
      if (oactive) then
         call choleski(  n, 
     &        dbl_mb(adrs), int_mb(adrma),  
     &        int_mb(adris), 
     &        dbl_mb(adrip),  info)  
C        if (info .ne. 0)
C    $        call ga_error('ga_chol: info non-zero ', info)
      endif
c     
      call ga_sync()            ! To reduce message buffering problems      
c
c***  store the Choleski factor matrix in global array g_s
c
      call ga_zero(g_s)
      if (oactive) then
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. myproc) then
              myelem = n-k+1
               call ga_put(g_s,k,n,k,k,dbl_mb(adrs+istart),1)
               istart = istart + myelem
            endif
         enddo
c     
c***     deallocate work arrays 
c     
         status = ma_pop_stack(hip) 
         status = ma_pop_stack(his) 
         status = ma_pop_stack(hs) 
         status = ma_pop_stack(hma) 
      endif
c
C#ifdef DATA_SERVER
C      FIX ME NOW
C#endif
c
c     Make sure that even processes that did not participate in the
c     process get info
c
c
      call ga_brdcst(Msg_ga_chol, info, mitob(1), 0)
      call ga_sync()
c
      end


