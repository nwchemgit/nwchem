      logical function spcart_init(lmaxin,normalize,all_spherical)
c $Id: int_spcart.F,v 1.2 1995-12-20 22:50:22 d3e129 Exp $
* 
* initialization of spherical to cartesean tranformation array
*... stored on heap.  
*... stored up to lmax values
*...
*   spcart(iccart,icsp,l)  => spcart((lmax+1)*(lmax+2)/2,1:2*lmax+1,0:lmax)
*    lmax = 5 h functions  => size = 21*11*6 = 1386 
*   store array 34% sparse for simplicity.  1386 doubles is 11 Kbytes.
*
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "spcartP.fh"
c::passed
      integer lmaxin    ! [input] init transformed values up to lmaxin
      logical normalize ! [input] normalize the coefficients for integral
*                                 transformations
*                         true for integral transformations
*                         false for ECP integral computations
      logical all_spherical ! [input] generate all spherical compoents eg 6 d sphericals etc.
*.
      external bd_spcart  ! needed for cray T3D to link properly
*
      logical spcart_terminate 
      external spcart_terminate 
c::local
      integer size_sp2c         ! size of array
      integer l_block_size      ! size of compressed array
      integer lval, l2, ls
c
*rak:: temporary
      if (all_spherical) call errquit
     &    ('spcart_init: all spherical components not working yet',911)
c
      if (sph_cart_init.eq.SPH_CART_INIT_VALUE) then
        if (lmaxin.gt.lmax_init) then
          if (.not.spcart_terminate()) call errquit
     &        ('spcart_init: error terminatine old spcart_init',911)
        else
          spcart_init = .true.
          return   ! initialization already done to cover lmaxin
        endif
      endif
*
      if (all_spherical) then
        sph_cart_allsph = .true.
      else
        sph_cart_allsph = .false.
      endif
c
      size_sp2c = lmaxin+1
      size_sp2c = size_sp2c*(2*lmaxin+1)
      size_sp2c = size_sp2c*(((lmaxin+1)*(lmaxin+2))/2)
c
      active_sp2c = 
     &    ma_alloc_get(mt_dbl,size_sp2c,'sph 2 cart trans array',
     &    h_sp2c,k_sp2c)
      if (.not.active_sp2c)  call errquit
     &    ('spcart_init: alloc_get failed for size',size_sp2c)
c
* generate transformation matricies by recursion
      call xlmcoeff(lmaxin,dbl_mb(k_sp2c),normalize)
* generate all spherical components
      call xlm_coeff_full(lmaxin,dbl_mb(k_sp2c),normalize)
* allocate memory for index array
      active_sp2c_lindx = ma_alloc_get(
     &    mt_int,(lmaxin+1),' ptrs array xlm sph 2 cart ',
     &    h_sp2c_lindx,k_sp2c_lindx)
      if (.not.active_sp2c_lindx) call errquit
     &    ('spcart_init: alloc_get failed (index) ',911)
* determine size of compressed transformation arrays
      l_block_size = 0
      do 00100 lval=0,Lmaxin
        l2 = (((lval+1)*(lval+2))/2)
        ls = (2*lval+1) 
        l_block_size = l_block_size + l2*ls
00100 continue
* allocate memory for compressed transformation arrays
      active_sp2c_cmp = ma_alloc_get
     &    (mt_dbl,l_block_size,'sph 2 cart trans array cmp',
     &    h_sp2c_cmp,k_sp2c_cmp)
      if (.not. active_sp2c_cmp) call errquit
     &    ('spcart_init: alloc_get failed (array cmp) ',911)

* set up pointers and copy recursion array to compressed transformation arrays
      call xlm_ptrs(lmaxin,dbl_mb(k_sp2c),dbl_mb(k_sp2c),
     &    dbl_mb(k_sp2c_cmp),l_block_size,int_mb(k_sp2c_lindx))
c
* free up recursion copy of transformation matricies
      if (.not.ma_free_heap(h_sp2c)) call errquit
     &    ('spcart_init: free heab failed (array) ',911)
      active_sp2c = .false.
      k_sp2c = 0
      h_sp2c = 0
c
      sph_cart_init = Sph_Cart_Init_Value
      lmax_init = lmaxin
      spcart_init = .true.
c
      end
c----------------------------------------------------------------------
      logical function spcart_terminate()
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"
c
c terminates spcart data structure and initialization
c
      if (sph_cart_init.eq.SPH_CART_INIT_VALUE) then
        spcart_terminate = .true.
        if (active_sp2c) then
          spcart_terminate = spcart_terminate .and.
     &        ma_free_heap(h_sp2c)
          active_sp2c = .false.
          k_sp2c = 0
          h_sp2c = 0
        endif
        if (active_sp2c_cmp) then
          spcart_terminate = spcart_terminate .and.
     &        ma_free_heap(h_sp2c_cmp)
          active_sp2c_cmp = .false.
          k_sp2c_cmp = 0
          h_sp2c_cmp = 0
        endif
        if (active_sp2c_lindx) then
          spcart_terminate = spcart_terminate .and.
     &        ma_free_heap(h_sp2c_lindx)
          active_sp2c_lindx = .false.
          k_sp2c_lindx = 0
          h_sp2c_lindx = 0
        endif
        if (.not.spcart_terminate) call errquit
     &      (' error freeing heap in spcart_terminate',911)
        sph_cart_init = 0
      else
        spcart_terminate = .false.
      endif
      end
c----------------------------------------------------------------------
      subroutine xlm_coeff_full(Ld,D,normalize)
      implicit none
*   not implemented yet 
      integer Ld
      double precision D(*)
      logical normalize
      end
c----------------------------------------------------------------------
* set up pointer information
      subroutine xlm_ptrs(Ld,D,Dall,Dcmp,ldcmp,Dindex)
#include "stdio.fh"
#include "mafdecls.fh"
#include "spcartP.fh"
c
c::passed
      integer Ld        ! [input] Lmax for how spcart_* was initialized
      integer ldcmp     ! [input] length of compressed transformation arrays
      double precision D((((Ld+1)*(Ld+2))/2),-Ld:Ld,0:Ld) ! [input] transformation matrix
      double precision Dall((((Ld+1)*(Ld+2))/2),
     &    (((Ld+1)*(Ld+2))/2),0:Ld) ! [input] transformation matrix
c
      double precision Dcmp(ldcmp)  ! [output] compressed transformation arrays
      integer Dindex(0:Ld)          ! [output] pointer for lth transform array in compressed array
c
c::local      
      integer lval
c
      if (sph_cart_allsph) call errquit
     &    ('xlm_ptrs: all spherical components not working yet',911)
c
c
      icount = 0
      do 00100 lval = 0,Ld
        Dindex(lval) = k_sp2c_cmp + icount  ! set pointer in index array
        do 00200 isp = -lval,lval
          do 00300 icart = 1,(((lval+1)*(lval+2))/2)
            icount = icount + 1
            Dcmp(icount) = D(icart,isp,lval)  ! form separate D(xyz,sph) arrays
00300     continue
00200   continue
00100 continue
      call spcart_print_both(D,ld)
      end
c----------------------------------------------------------------------
      subroutine spcart_print_both(D,ld)
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"
      integer ld
      double precision D((((Ld+1)*(Ld+2))/2),-Ld:Ld,0:Ld) ! [input] transformation matrix
c
      integer lval,l2s, ic, is
      double precision diff
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
      do lval = 0,ld
        write(6,*)' d matrix '
        l2s = (lval+1)*(lval+2)/2
        do is = -lval,lval
          do ic = 1,l2s
            diff = D(ic,is,lval)-Dtrans(ic,is,lval)
            write(6,*) lval,is,ic,D(ic,is,lval),Dtrans(ic,is,lval),diff
          enddo
        enddo
      enddo
      end
c----------------------------------------------------------------------
      integer function spcart_trns_ptr(lval,lcart,lsp)
c return pointer in ma to transformation matrix for lval
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"

      integer lval  ! [input] l-value of requested matrix
      integer lcart ! [output] number of cartesian components
      integer lsp   ! [output] number of spherical components
c
      spcart_trns_ptr = int_mb(k_sp2c_lindx+lval)
      lcart = ((lval+1)*(lval+2))/2
      if (sph_cart_allsph) then
        lsp = lcart
      else
        lsp = 2*lval+1
      endif
      end
c----------------------------------------------------------------------
      Block data bd_spcart

#include "spcartP.fh"

      data sph_cart_init /0/
      data lmax_init     /0/
      data h_sp2c        /0/
      data k_sp2c        /0/
      data h_sp2c_cmp    /0/
      data k_sp2c_cmp    /0/
      data h_sp2c_lindx  /0/
      data k_sp2c_lindx  /0/
      data sph_cart_allsph   /.false./
      data active_sp2c       /.false./
      data active_sp2c_cmp   /.false./
      data active_sp2c_lindx /.false./
      end
c----------------------------------------------------------------------
      subroutine spcart_a_s(blockin, blockout, ndima, ls, in_place)
      implicit none
c
c  transforms a block of integrals with the Ls function is the slowest dimension
c  e.g., blockin(ndima,L2s)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ls     ! [input] angular momentum of block
      double precision blockin (ndima,((ls+1)*(ls+2)/2))  ! [input] matrix 
*      double precision blockout(ndima,(2*ls+1))           ! [output] matrix 
      double precision blockout(ndima,-ls:ls)           ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i,j,k
      integer L2s
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (ls.lt.2) then
c...        ((ls+1)*(ls+2)/2)) = 2*ls + 1
        call dcopy((ndima*(2*ls+1)),blockin,1,blockout,1)
        return
      endif
c

      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
#if defined(USE_DGEMM)
c
c dgemm call for blockin(ndima,l2s)*dtrans(l2s,2l+1) = blockout(ndima,2l+1)
c
      L2s = ((ls+1)*(ls+2)/2)
      call dgemm('n','n',ndima,(2*ls+1),L2s,1.0d00,blockin,ndima,
     &    dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,0.0d00,blockout,ndima)
#else
c
c  blockout(ndima,2l+1) = blockin(ndima,l2s)*d_spcart(l2s,2l+1,lp=ls)
c      
      L2s = ((ls+1)*(ls+2)/2)
      do 00100 i=1,ndima
        do 00200 j=(-ls),ls
          blockout(i,j) = 0.0d00
          do 00300 k = 1,L2s
            blockout(i,j) = blockout(i,j) + blockin(i,k)*Dtrans(k,j,Ls)
00300     continue
00200   continue
00100 continue
#endif
      end
c----------------------------------------------------------------------
      subroutine spcart_s_a(blockin, blockout, ndima, ls, in_place)
      implicit none
c
c  transforms a block of "integrals" with the ls function is the leading dimension
c  e.g., blockin(L2s,ndima)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ls     ! [input] angular momentum of block
      double precision blockin (((ls+1)*(ls+2)/2),ndima)  ! [input] matrix 
*      double precision blockout((2*ls+1),ndima)           ! [output] matrix 
      double precision blockout(-ls:ls,ndima)           ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i,j,k
      integer L2s
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      write(6,*)'s_a,(l2s,ndima) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (ls.lt.2) then
c...        ((ls+1)*(ls+2)/2)) = 2*ls + 1
        call dcopy((ndima*(2*ls+1)),blockin,1,blockout,1)
        return
      endif
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
#if defined(USE_DGEMM)
c
c dgemm call for:
c  Transpose(d_spcart(l2s,2l+1))*blockin(l2s,ndima) = blockout(2l+1,ndima)
c
      L2s = ((ls+1)*(ls+2)/2)
      call dgemm('t','n',(2*ls+1),ndima,L2s,1.0d00,
     &      dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,
     &      blockin,L2s,0.0d00,blockout,(2*ls+1))
#else
c
c blockout(2l+1,ndima) = blockin(l2s,ndima)*d_spcart(l2s,2l+1,lp=ls)
c      
      L2s = ((ls+1)*(ls+2)/2)
      call dfill(((2*ls+1)*ndima),0.0d00,blockout,1)
      do 00100 i=(-ls),ls
        do 00200 j=1,ndima
          do 00300 k = 1,L2s
*rak:            write(7,*)i, j, k, blockin(k,j), Dtrans(k,i,Ls),
*rak:     &          blockout(i,j)
            blockout(i,j) = blockout(i,j) + blockin(k,j)*Dtrans(k,i,Ls)
*rak:            write(7,*)blockout(i,j)
00300     continue
00200   continue
00100 continue
#endif
      end
c----------------------------------------------------------------------
      subroutine spcart_a_s_b(blockin, blockout, ndima, ndimb, ls,
     &    in_place)
      implicit none
c
c  transforms a block of "integrals" with the ls function is ordered between 
c  a leading dimension and trailing dimension.
c  e.g., blockin(nidima,L2s,ndimb)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] tailing dimension of block
      integer ls     ! [input] angular momentum of block
      double precision blockin (ndima,((ls+1)*(ls+2)/2),ndimb) ! [input] matrix 
*      double precision blockout(ndima,(2*ls+1),ndimb)          ! [output] matrix 
      double precision blockout(ndima,-ls:ls,ndimb)          ! [output] matrix 
      logical in_place  ! [input] true if blockin and blockout are the same pointer
c: local
      integer i,j,k, m
      integer L2s
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      write(6,*)'a_s,(ndima,l2s,ndimb) ndima = ',ndima,
*rak:     &      '   ls = ',ls,'  ndimb = ',ndimb
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
      if (ls.lt.2) then
c...        ((ls+1)*(ls+2)/2)) = 2*ls + 1
        call dcopy((ndima*ndimb*(2*ls+1)),blockin,1,blockout,1)
        return
      endif
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
#if defined(USE_DGEMM)
c
c dgemm for all ndimb
c
      L2s = ((ls+1)*(ls+2)/2)
      do 00100 m = 1,ndimb
        call dgemm('n','n',ndima,(2*ls+1),L2s,
     &        1.0d00,blockin(1,1,m),ndima,
     &        dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,
     &        0.0d00,blockout(1,1,m),ndima)
00100 continue
#else
      L2s = ((ls+1)*(ls+2)/2)
      do 00100 m = 1,ndimb
        do 00200 i=1,ndima
          do 00300 j=(-ls),ls
            blockout(i,j,m) = 0.0d00
            do 00400 k = 1,L2s
              blockout(i,j,m) = blockout(i,j,m) +
     &            blockin(i,k,m)*Dtrans(k,j,Ls)
00400       continue
00300     continue
00200   continue
00100 continue
#endif        
      end
