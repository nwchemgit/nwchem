      subroutine edo_diag_compl(g_ar, g_ac, g_vr, g_vc, eval)
c     $Id: ga_compl.F,v 1.3 2000-10-25 22:50:14 windus Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer g_ar,g_ac         ! matrix to diagonalize (real and imag)
      integer g_vr,g_vc         ! global matrix to return evecs
      double precision eval(*)  ! local array to return evals
      double precision ee,dotres,toll, toll1
      parameter(toll=1.d-9,toll1=0.99)
c     
c     
c     solve the complex symmetric eigenvalue problem 
c     for a hermitian complex matrix mapping it to
c     a real one in the following way
c     ( AR  -AC)  (VR VR)  = (VR VR) (E  )
c     ( AC   AR)  (VC VC)  = (VC VC) (  E)
c     returning
c     all eigenvectors and values in ascending order
c     
c     the input matrix is not destroyed
c     
c     proc - # of processors
c     pans - panel size in columns
c     n    - problem size
c     
c     
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     
      integer iproc, pan, istart, istart1, iend, j,k,k1,n, nn, 
     *     npan, proc, myproc
      integer mypanel, myelem, mycol, elemz, ipan, ik, info
      logical status
      logical loop1,loop2
c     
      integer hma, adrma        !mapa
      integer hmz, adrmz        !mapz
      integer ga_diag_type,sndvec,brdres
      integer lenmes
c     
      integer ha, adra          !a
      integer hz, adrz          !z
c     
      integer isize, rsize, ptr_size
      integer his, adris
      integer hsc, adrsc
      integer hip, adrip
      integer hj,adrj
      integer he,adre
      integer node1,node2
      integer igua,imazz
c     
      external fmemreq, pdspev
c     
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimV1, dimV2, typeV
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
      call ga_check_handle(g_ar, 'ga_diag: ar')
      call ga_check_handle(g_vr, 'ga_diag: vr')
      call ga_check_handle(g_ac, 'ga_diag: ac')
      call ga_check_handle(g_vc, 'ga_diag: vc')
      call ga_inquire(g_ar, typeA, dimA1, dimA2)
      call ga_inquire(g_vr, typeV, dimV1, dimV2)
c     
c     
      if (dimA1 .ne. dimA2) then
         call ga_error('ga_diag_compl: can only diag square matrix ', 1)
      else if (dimV1 .ne. dimV2) then
         call ga_error('ga_diag_compl: can only diag square matrix ', 3)
      else if(dimA1.ne.dimV1 )then
         call ga_error('ga_diag_compl: A & V dimension mismatch ', 4)
      endif
c     
c     call mxpara(ienv)
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
c     print *,myproc,' participate in eig: ',proc
c     call ffflush(6)
      oactive = myproc .lt. proc
      if (oactive) then
         pan = 1                !panel size
         npan = n/pan
c     
c     allocate map arrays
c     
         status = ma_push_get(mt_int,n,'mapa',hma,adrma)
         status = ma_push_get(mt_int,n,'mapz',hmz,adrmz)
c     
c     determine distribution of a and z
c     
         mypanel = 0            ! number of panels at given processor
         mycol   = 0            ! number of columns at given processor
         myelem  = 0            ! number of elements at given processor
c     
c     allocate first half of the panels - forwards
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
c     allocate second half of the panels - backwards
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
c***  actually, there is one more panel left for the mismatch 
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
c***  allocate a and z arrays
c     
         status = ma_push_get(mt_dbl,myelem,'a',ha,adra)
         status = ma_push_get(mt_dbl,elemz,'z',hz,adrz)
c     
c***  load the fock (A) matrix
c     [ AR    ]
c     [ AC  AR] 
c     
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. myproc) then
               myelem = n-k+1
               if(k.le.nn) then
                  call ga_get(
     *                 g_ar,k,nn,k,k,dbl_mb(adra+istart),1)
                  call ga_get(
     *                 g_ac,1,nn,k,k,dbl_mb(adra+istart+(nn-k+1)),1)
               else
                  k1=k-nn
                  call ga_get(
     *                 g_ar,k1,nn,k1,k1,dbl_mb(adra+istart),1)
               endif
               istart = istart + myelem
            endif
         enddo
c     
c***  allocate some memory for compl eigensolver 
c     
         status = ma_push_get(mt_int,n,'iscratch',his,adris)
         call fmemreq(1,n,int_mb(adrma),int_mb(adrma),int_mb(adrmz),
     &        isize, rsize, ptr_size,int_mb(adris))
         status = ma_pop_stack(his)
         if(.not. ma_push_get(mt_int,isize,'iscratch',his,adris))
     &        call ga_error('ga_diag_compl allocator',1)
         if(.not. ma_push_get(mt_dbl,rsize,'scratch',hsc,adrsc))
     &        call ga_error('ga_diag_compl allocator',2)
         if(.not. ma_push_get(mt_dbl,ptr_size,'iptr',hip,adrip))
     &        call ga_error('ga_diag_compl allocator',3)
c     
      endif
c     
c     
      call ga_sync()            ! To reduce message buffering problems
c     
c***  calling PEIGS parallel eigensolver
c     
      if (oactive) then
         call pdspev(  n, dbl_mb(adra), int_mb(adrma),
     &        dbl_mb(adrz), int_mb(adrmz),  
     &        eval,  
     &        int_mb(adris), isize,  
     &        dbl_mb(adrip), ptr_size,  
     &        dbl_mb(adrsc), rsize, info)  
         if (info .ne. 0)
     $        call ga_error('ga_diag_compl: info non-zero ', info)
      endif
c     
      call ga_sync()            ! To reduce message buffering problems

      if (oactive) then
c     
c     discard redundant eigenvals and eigenvecs
c     
         status = ma_push_get(mt_int,n,'junk',hj,adrj)
         status = ma_push_get(mt_dbl,n,'junke',he,adre)
         do k = 1,44
            write(*,'(11(1x,f6.3))')(dbl_mb(adrz+44*(k-1)+k1),k1=0,43)
            write(*,*)
         enddo 
         do k=1,n
            int_mb(adrj+k-1)=0
c     write(0,*) ' k eval ',k,eval(k)
         enddo
         istart=0
         imazz = 2              ! take care of compiler warnings
         do k = 1,n
            loop1=.false.
            if(int_mb(-1+adrmz+k).eq.myproc) loop1=.true.
            if (int_mb(adrj+k-1).eq.0) then
               ee=eval(k)
               igua=0
               do k1=k+1,n
                  if(abs(eval(k1)-ee).lt.toll) then
                     igua=igua+1
                     imazz=k1
                  endif
               enddo
               if(igua.eq.1) then
                  int_mb(adrj+imazz-1)=1
                  goto 123
               endif
               istart1=istart
               if(int_mb(-1+adrmz+k).eq.myproc) istart1=istart1+n
               do k1 = k+1,n
                  loop2=.false.
                  if(int_mb(-1+adrmz+k1).eq.myproc) loop2=.true.
                  if(k.ne.k1.and.abs(eval(k1)-ee).lt.toll) then
c     
c     dot product
c     real + imaginary part
c     
                     node1=int_mb(adrmz+k-1)
                     node2=int_mb(adrmz+k1-1)
c     
c     send everything to node1
c     
                     sndvec=65536+137
                     if(node1.ne.node2) then
                        if(myproc.eq.node2) 
     *                       call snd(sndvec,dbl_mb(adrz+istart1),
     *                       mdtob(n),node1,1)
                        if(myproc.eq.node1) 
     *                       call rcv(sndvec,dbl_mb(adre),mdtob(n),
     *                       lenmes,-1,node2,1)
                     else
                        if(myproc.eq.node1) then
                           do j=0,n-1
                              dbl_mb(adre+j)=dbl_mb(adrz+istart1+j)
                           enddo
                        endif
                     endif
                     if(myproc.eq.node1) then
                        dotres=0.d0
                        do j=0,nn-1
                           dotres=dotres+
     *                          dbl_mb(adrz+istart+j+nn)*dbl_mb(adre+j)+
     *                          dbl_mb(adrz+istart+j)*dbl_mb(adre+j+nn)
                        enddo
                     endif
                     brdres=65536+23+k+k1
                     if (myproc.eq.node1) then
                        do j=0,proc-1
                           if(j.ne.myproc)
     *                          call snd(brdres,dotres,mdtob(1),j,1)
                        enddo
                     else
                        call rcv(brdres,dotres,mdtob(1),
     *                       lenmes,-1,node1,1)
                     endif
                     if (abs(dotres).gt.toll) then
                        int_mb(adrj+k1-1)=1
                        goto 123
                     endif
                  endif
                  if(loop2) istart1=istart1+n
               enddo
 123           continue
            endif
            if(loop1) istart=istart+n
         enddo
c     
c     discard unuseful eigenvec
c     
         do k=1,n
            dbl_mb(adre+k-1)=eval(k)
         enddo
         k1=0
         istart = 0
         do k=1,n
            loop1=.false.
            if( int_mb(-1+adrmz+k) .eq. myproc) loop1=.true.
            if(int_mb(adrj+k-1).eq.0) then
               k1=k1+1
               if(k1.gt.nn) then
                  if (myproc.eq.0) then
                     print *,' K1 ',k1
                     write(0,*) ' k evals IJUNK'
                     do j=1,n
                        write(0,*) j,dbl_mb(adre-1+j),int_mb(adrj+j-1)
                     enddo
                  endif
                  call ga_sync
                  call ga_error('ga_diag_compl: 
     *                 failed to reduce eigenvec', 1)
               endif
               eval(k1)=dbl_mb(adre+k-1)
c     
c***  store the eigenvector matrix in global array g_vr and g_vc
c     
               if(loop1) then
                  call ga_put(g_vr,1,nn,k1,k1,dbl_mb(adrz+istart),1)
                  call ga_put(g_vc,1,nn,k1,k1,dbl_mb(adrz+istart+nn),1)
               endif
            endif
            if(loop1) istart = istart + n
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
         status = ma_pop_stack(hmz) 
         status = ma_pop_stack(hma) 
      endif
c     
c     Make sure that even processes that did not participate in the
c     diagonalization have the eigen values
c     
      ga_diag_type = 11
      call brdcst(ga_diag_type, eval, mdtob(n), 0)
c     
#ifdef DATA_SERVER
      FIX ME NOW
#endif
c     
      call ga_sync()
c     
      end
