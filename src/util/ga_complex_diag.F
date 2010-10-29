      subroutine ga_diag_compl(g_ar, g_ac, g_vr, g_vc, eval)
*
* $Id$
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
      integer g_ar,g_ac         ! matrix to diagonalize (real and imag)
      integer g_vr,g_vc         ! global matrix to return evecs
      double precision eval(*)  ! local array to return evals
      double precision ee,dotres,toll
      parameter(toll=1.d-9)
      intrinsic abs
c     
c       
c     solve thecomplex symmetric eigenvalue problem 
c     for a hermitiancomplex matrix mapping it to
c     a real one in the following way
c     
c     
c     
c     ( AR  -AC)  (VR VR)  = (VR VR) (E  )
c     ( AC   AR)  (VC VC)  = (VC VC) (  E)
c
c since AR + i AC is hermitian this is a symmetric representation
c
c       returning
c       all eigenvectors and values in ascending order
c       
c       the input matrix is not destroyed
c       
c       proc - # of processors
c       pans - panel size incolumns
c       n    - problem size
c       
c       
c       
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       
	integer iproc, pan, istart, istart1, iend, j,k,k1,n, nn, 
     $     npan, proc, myproc
	integer mypanel, myelem, mycol, elemz, ipan, ik, info
	logical status
	logical loop1,loop2
c       
	integer hma, adrma	!mapa
	integer hmz, adrmz	!mapz
	integer ga_diag_type,sndvec,brdres
	integer lenmes
c
	integer ha, adra	!a
	integer hz, adrz	!z
c       
	integer isize, rsize, ptr_size
	integer his, adris
	integer hsc, adrsc
	integer hip, adrip
	integer hj,adrj
	integer he,adre
	integer heval,adreval
	integer node1,node2
	integer igua,imazz
c       
	external fmemreq, pdspev
c       
	logical oactive		! true iff this process participates
	character*80 name
	integer dimA1, dimA2, typeA
	integer dimV1, dimV2, typeV
c       
	integer Elem
	Elem(istart,iend,n)=((iend-istart+1)*(2*n-istart-iend+2))/2
c       
c       
c***    check environment
c       
	myproc = ga_nodeid()
	proc   = ga_nnodes()
c       
c***    check GA info
c       
	call ga_check_handle(g_ar, 'ga_diag: ar')
	call ga_check_handle(g_vr, 'ga_diag: vr')
	call ga_check_handle(g_ac, 'ga_diag: ac')
	call ga_check_handle(g_vc, 'ga_diag: vc')
	call ga_inquire(g_ar, typeA, dimA1, dimA2)
	call ga_inquire(g_vr, typeV, dimV1, dimV2)
c       
c       
	if (dimA1 .ne. dimA2) then
	   call ga_error('ga_diag_compl:can only diag square matrix ',2
     $          )
	else if (dimV1 .ne. dimV2) then
	   call ga_error('ga_diag_compl:can only diag square matrix ',3
     $          )
	else if(dimA1.ne.dimV1 )then
         call ga_error('ga_diag_compl: A & V dimension mismatch ', 4)
	endif
c       
c       call mxpara(ienv)
	nn = dimA1
	n = 2*nn
c       
	call ga_sync()
c
       if (n .lt. 350) then
          proc = min(8,ga_nnodes())
       else if (n .lt. 700) then
          proc = min(16,ga_nnodes())
       else if (n .lt. 1000) then
          proc = min(20,ga_nnodes())
       else if (n .lt. 2000) then
          proc = min(128,ga_nnodes())
       else
          proc = min(256,ga_nnodes())
       endif
c
c         print *,myproc,' participate in eig: ',proc
c        call ffflush(6)
       oactive = myproc .lt. proc
       if (oactive) then
          pan = 1                !panel size
          npan = n/pan
c     
c        allocate map arrays
c
          status = ma_push_get(mt_int,n,'mapa',hma,adrma)
          status = ma_push_get(mt_int,n,'mapz',hmz,adrmz)
c         
c        determine distribution of a and z
c     
          mypanel = 0            ! number of panels at given processor
          mycol   = 0            ! number ofcolumns at given processor
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
            int_mb(-1+adrmz+k) = iproc
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
                int_mb(-1+adrmz+k) = iproc
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
             int_mb(-1+adrmz+k) = iproc
          enddo
          mycol = mypanel * pan
          if(iproc .eq. myproc)then
             myelem = myelem + elem(istart,iend,n)
             mycol  = mycol + iend - istart + 1	
          endif
          elemz = mycol * n
c
c***     allocate a and z arrays
c     
          status = ma_push_get(mt_dbl,n,'heval',heval,adreval)
          if(.not. status)
     $         call ga_error('ga_diag_compl allocator heval',1)
c     
          status = ma_push_get(mt_dbl,myelem,'a',ha,adra)
          status = ma_push_get(mt_dbl,elemz,'z',hz,adrz)
c
c***     load the fock (A) matrix
c        [ AR    ]
c        [ AC  AR] 
c
          istart = 0
          do k = 1, n
             if( int_mb(-1+adrma+k) .eq. myproc) then
                myelem = n-k+1
               if(k.le.nn) then
                call ga_get(
     $                 g_ar,k,nn,k,k,dbl_mb(adra+istart),1)
                call ga_get(
     $               g_ac,1,nn,k,k,dbl_mb(adra+istart+(nn-k+1)),1)
               else
                 k1=k-nn
                call ga_get(
     $                g_ar,k1,nn,k1,k1,dbl_mb(adra+istart),1)
               endif
               istart = istart + myelem
             endif
          enddo
c     
c***     allocate some memory forcompl eigensolver 
c
          status = ma_push_get(mt_int,n,'iscratch',his,adris)
          call fmemreq(1,n,int_mb(adrma),int_mb(adrma),int_mb(adrmz),
     $         isize, rsize, ptr_size,int_mb(adris))
          status = ma_pop_stack(his)
          if(.not. ma_push_get(mt_int,isize,'iscratch',his,adris))
     $         call ga_error('ga_diag_compl allocator',1)
          if(.not. ma_push_get(mt_dbl,rsize,'scratch',hsc,adrsc))
     $         call ga_error('ga_diag_compl allocator',2)
          if(.not. ma_push_get(mt_dbl,ptr_size,'iptr',hip,adrip))
     $         call ga_error('ga_diag_compl allocator',3)
c     
       endif
c
c
      call ga_sync()            ! To reduce message buffering problems
c
c***  calling PEIGS parallel eigensolver
c     
      if (oactive) then
c
c set info to use complex diag
c
         info = -1
         call pdspev(  n, dbl_mb(adra), int_mb(adrma),
     $        dbl_mb(adrz), int_mb(adrmz),  
     $        dbl_mb(adreval),  
     $        int_mb(adris), isize,  
     $        dbl_mb(adrip), ptr_size,  
     $        dbl_mb(adrsc), rsize, info)
         if (info .ne. 0)
     $        call ga_error('ga_diag_compl: info non-zero ', info)
         call pdcomplex(n, dbl_mb(adrz), int_mb(adrmz),
     $        dbl_mb(adreval),
     $        dbl_mb(adrsc), dbl_mb(adris), info)
         if (info .ne. 0)
     $        call ga_error(
     $        'ga_diag_comp: pdcomplex info non-zero ', info)
      endif
c     
c     to reduce out of order message problems between ga and non-ga process
c     
c     
      call ga_sync
c     
c***  store the eigenvector matrix in global array g_vr and g_vc
c     
      n = nn
      if (oactive) then
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrmz+k) .eq. myproc) then
               call ga_put(g_vr,1,nn,k,k,dbl_mb(adrz+istart),1)
               call ga_put(g_vc,1,nn,k,k,dbl_mb(adrz+istart+nn),1)
               istart = istart + n
            endif
      enddo
c     
c     set n = n/2
c
      n = nn
c     
c***  deallocate work arrays 
c     
      status = ma_pop_stack(he) 
      status = ma_pop_stack(hj) 
      status = ma_pop_stack(hip) 
      status = ma_pop_stack(hsc) 
      status = ma_pop_stack(his) 
      status = ma_pop_stack(hz) 
      status = ma_pop_stack(ha) 
      status = ma_pop_stack(heval)
      status = ma_pop_stack(hmz) 
      status = ma_pop_stack(hma) 
      endif

c     
c     Make sure that even processes that did not participate in the
c     diagonalization have the eigen values
c     
c     
#ifdef DATA_SERVER
      FIX ME NOW
#endif
c     
      call ga_sync()
c     
      end
      
