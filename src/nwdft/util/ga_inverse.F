      subroutine ga_inverse(g_a, g_s)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer g_a               ! input matrix
      integer g_s               ! output matrix
c
c
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
c
      integer hma, adrma        !mapa
c
      integer hs, adrs          !s
c
      integer his, adris
      integer hip, adrip
c
      external  inverseL,fmemreq
c
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimS1, dimS2, typeS
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
      call ga_check_handle(g_a, 'ga_inverse: a')
      call ga_check_handle(g_s, 'ga_inverse: s')
      call ga_inquire(g_a,  typeA, dimA1, dimA2)
      call ga_inquire(g_s,  typeS, dimS1, dimS2)
c
c
      if (dimA1 .ne. dimA2) then
         call ga_error('ga_inverse: can only invert square matrix ', 1)
      else if (dimS1 .ne. dimS2) then
         call ga_error('ga_inverse: can only invert square matrix ', 3)
      else if(dimA1.ne.dimS1 )then
         call ga_error('ga_inverse: A & S dimension mismatch ', 4)
      endif
c
      n = dimA1
      call ga_sync()
      proc = min(n/15, ga_nnodes())
c     if n.lt.15 proc becomes zero
      proc=max(1,proc)
      oactive = myproc .lt. proc
      if (oactive) then
         pan = 1                !panel size
         npan = n/pan
c     
c        allocate map arrays
c
         if(.not.ma_push_get(mt_int,n,'mapa',hma,adrma))
     &     call errquit('ga_inverse: cannot allocate mapa',0, MA_ERR)
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
         if(.not.ma_push_get(mt_dbl,elems,'s',hs,adrs))
     &     call errquit('ga_inverse: cannot allocate s',0, MA_ERR)
c
c***     load the fock (A) matrix to (S)
c
         
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. myproc) then
               myelem = n-k+1
               call ga_get(g_a,k,n,k,k,dbl_mb(adrs+istart),n-k+1)
               istart = istart + myelem
            endif
         enddo
c     
c***     allocate some memory 
c

         if(.not.ma_push_get(mt_int,(6*n),'iwork',his,adris))
     &     call errquit('ga_inverse: cannot allocate iwork',0, MA_ERR)

         if(.not.ma_push_get(mt_dbl,(6*n+1),'work',hip,adrip))
     &     call errquit('ga_inverse: cannot allocate work',0, MA_ERR)
     
      endif
c
c
      call ga_sync()            ! To reduce message buffering problems
c
c***  calling PEIGS parallel inverse
c
      if (oactive) then
         call inverseL(  n, 
     &        int_mb(adrma),  dbl_mb(adrs), 
     &        int_mb(adris), 
     &        dbl_mb(adrip),  info)  
         if (info .ne. 0)
     $        call ga_error('ga_inverse: info non-zero ', info)
      endif
c     
      call ga_sync()            ! To reduce message buffering problems      
c
c***  store the Inverse matrix in global array g_s
c
      call ga_zero(g_s)
      if (oactive) then
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. myproc) then
              myelem = n-k+1
               call ga_put(g_s,k,n,k,k,dbl_mb(adrs+istart),n-k+1)
               istart = istart + myelem
            endif
         enddo
c     
c***     deallocate work arrays 
c     
         if(.not.ma_pop_stack(hip))
     &     call errquit('ga_inverse: cannot pop stack',0, MA_ERR)

         if(.not.ma_pop_stack(his))
     &     call errquit('ga_inverse: cannot pop stack',0, MA_ERR)

         if(.not.ma_pop_stack(hs))
     &     call errquit('ga_inverse: cannot pop stack',0, MA_ERR)

         if(.not.ma_pop_stack(hma))
     &     call errquit('ga_inverse: cannot pop stack',0, MA_ERR)
      endif
c
      call ga_sync()
c
      end


