      logical function spcart_init(lmaxin,normalize,all_spherical)
c $Id: int_spcart.F,v 1.10 1996-03-20 00:53:50 d3e129 Exp $
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
*.......................................................................
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
*.......................................................................
      subroutine xlm_coeff_full(Ld,D,normalize)
      implicit none
*   not implemented yet 
      integer Ld
      double precision D(*)
      logical normalize
      end
*.......................................................................
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
*     call spcart_print_both(D,ld)
      end
*.......................................................................
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
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_print_both: spcart not initialized properly',
     &    sph_cart_init)
c
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
*.......................................................................
      subroutine spcart_print_dtrans(ld)
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"
      integer ld
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_print_dtrans: spcart not initialized properly',
     &    sph_cart_init)
c
      write(6,*)' trans matrix used for Lval =',ld
      call output(dbl_mb((int_mb(k_sp2c_lindx+ld))),1,
     &    ((ld+1)*(ld+2)/2),1,(2*ld+1),
     &    ((ld+1)*(ld+2)/2),(2*ld+1),1)
      end
*.......................................................................
      integer function spcart_trns_ptr(lval,lcart,lsp)
c return pointer in ma to transformation matrix for lval
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"

      integer lval  ! [input] l-value of requested matrix
      integer lcart ! [output] number of cartesian components
      integer lsp   ! [output] number of spherical components
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_trns_ptr: spcart not initialized properly',
     &    sph_cart_init)
c
      spcart_trns_ptr = int_mb(k_sp2c_lindx+lval)
      lcart = ((lval+1)*(lval+2))/2
      if (sph_cart_allsph) then
        lsp = lcart
      else
        lsp = 2*lval+1
      endif
      end
*.......................................................................
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
*.......................................................................
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
      double precision blockout(ndima,-ls:ls)           ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
#define SPECIAL
#if defined(USE_LOOPS)
      integer i,j,k
#endif
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
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_s: spcart not initialized properly',
     &    sph_cart_init)
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
c
*      call spcart_print_Dtrans(ls)
c
#if defined(SPECIAL)
      if (ls.eq.2) then
        call spcart_a_sd(blockin,blockout,ndima,in_place)
        return
      else if (ls.eq.3) then
        call spcart_a_sf(blockin,blockout,ndima,in_place)
        return
      else if (ls.eq.4) then
        call spcart_a_sg(blockin,blockout,ndima,in_place)
        return
      endif
#endif
#if defined(USE_LOOPS)
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
#else
c
c
c dgemm call for blockin(ndima,l2s)*dtrans(l2s,2l+1) = blockout(ndima,2l+1)
c
      L2s = ((ls+1)*(ls+2)/2)
      call dgemm('n','n',ndima,(2*ls+1),L2s,1.0d00,blockin,ndima,
     &    dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,0.0d00,blockout,ndima)
c
#endif
      end
*.......................................................................
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
      double precision blockout(-ls:ls,ndima)           ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
#if defined(USE_LOOPS)
      integer i,j,k
#endif
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
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_s_a: spcart not initialized properly',
     &    sph_cart_init)
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
c
*      call spcart_print_Dtrans(ls)
c
#if defined(SPECIAL)
      if (ls.eq.2) then
        call spcart_sd_a(blockin,blockout,ndima,in_place)
        return
      else if (ls.eq.3) then
        call spcart_sf_a(blockin,blockout,ndima,in_place)
        return
      else if (ls.eq.4) then
        call spcart_sg_a(blockin,blockout,ndima,in_place)
        return
      endif
#endif
#if defined(USE_LOOPS)
c
c blockout(2l+1,ndima) = blockin(l2s,ndima)*d_spcart(l2s,2l+1,lp=ls)
c      
      L2s = ((ls+1)*(ls+2)/2)
      call dfill(((2*ls+1)*ndima),0.0d00,blockout,1)
      do 00100 i=(-ls),ls
        do 00200 j=1,ndima
          do 00300 k = 1,L2s
            blockout(i,j) = blockout(i,j) + blockin(k,j)*Dtrans(k,i,Ls)
00300     continue
00200   continue
00100 continue
#else
c
c dgemm call for:
c  Transpose(d_spcart(l2s,2l+1))*blockin(l2s,ndima) = blockout(2l+1,ndima)
c
      L2s = ((ls+1)*(ls+2)/2)
      call dgemm('t','n',(2*ls+1),ndima,L2s,1.0d00,
     &      dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,
     &      blockin,L2s,0.0d00,blockout,(2*ls+1))
#endif
      end
*.......................................................................
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
      double precision blockout(ndima,-ls:ls,ndimb)            ! [output] matrix 
      logical in_place  ! [input] true if blockin and blockout are the same pointer
c: local
#if defined(USE_LOOPS)
      integer i,j,k
#endif
      integer m
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
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_s_b: spcart not initialized properly',
     &    sph_cart_init)
c
      if (ls.lt.2) then
c...        ((ls+1)*(ls+2)/2)) = 2*ls + 1
        call dcopy((ndima*ndimb*(2*ls+1)),blockin,1,blockout,1)
        return
      endif
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
*      call spcart_print_Dtrans(ls)
c
#if defined(SPECIAL)
      if (ls.eq.2) then
        call spcart_a_sd_b(blockin,blockout,ndima,ndimb,in_place)
        return
      else if (ls.eq.3) then
        call spcart_a_sf_b(blockin,blockout,ndima,ndimb,in_place)
        return
      else if (ls.eq.4) then
        call spcart_a_sg_b(blockin,blockout,ndima,ndimb,in_place)
        return
      endif
#endif
c
#if defined(USE_LOOPS)
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
#else
c dgemm for all ndimb
c
      L2s = ((ls+1)*(ls+2)/2)
      do 00100 m = 1,ndimb
        call dgemm('n','n',ndima,(2*ls+1),L2s,
     &        1.0d00,blockin(1,1,m),ndima,
     &        dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,
     &        0.0d00,blockout(1,1,m),ndima)
00100 continue
#endif        
      end
*.......................................................................
      subroutine spcart_a_sd(blockin, blockout, ndima, in_place)
      implicit none
c
c  transforms a block of integrals where the d function is the slowest dimension
c  e.g., blockin(ndima,6) and blockout(ndima,5)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      double precision blockin (ndima,6)    ! [input] matrix 
      double precision blockout(ndima,-2:2) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i
*rak:      integer ls     ! [fixed at 2] angular momentum of block
      double precision dt2m2, dt5m1, dt10, dt40, dt60, dt31, dt12, dt42
      double precision bin1, bin2, bin3, bin4, bin5, bin6
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=2
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_sd: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c  spint(-2) = cint(2)*dtrans(2,-2,2)
c  spint(-1) = cint(5)*dtrans(5,-1,2)
c  spint( 0) = cint(1)*dtrans(1,0,2) + cint(4)*dtrans(4,0,2) + cint(6)*dtrans(6,0,2)
c  spint( 1) = cint(3)*dtrans(3,1,2)
c  spint( 2) = cint(1)*dtrans(1,2,2) + cint(4)*dtrans(4,2,2)
      dt2m2 = dtrans(2,-2,2)
      dt5m1 = dtrans(5,-1,2)
      dt10  = dtrans(1, 0,2)
      dt40  = dtrans(4, 0,2)
      dt60  = dtrans(6, 0,2)
      dt31  = dtrans(3, 1,2)
      dt12  = dtrans(1, 2,2)
      dt42  = dtrans(4, 2,2)
      do 00100 i=1,ndima
        bin1 = blockin(i,1)
        bin2 = blockin(i,2)
        bin3 = blockin(i,3)
        bin4 = blockin(i,4)
        bin5 = blockin(i,5)
        bin6 = blockin(i,6)
        blockout(i,-2) = bin2*dt2m2
        blockout(i,-1) = bin5*dt5m1
        blockout(i, 0) = bin1*dt10 +
     &                   bin4*dt40 +
     &                   bin6*dt60
        blockout(i, 1) = bin3*dt31
        blockout(i, 2) = bin1*dt12 +
     &                   bin4*dt42
00100 continue
      end
*.......................................................................
      subroutine spcart_sd_a(blockin, blockout, ndima, in_place)
      implicit none
c
c  transforms a block of integrals where the d function is the fastest dimension
c  e.g., blockin(6,ndima) and blockout(5,ndima)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      double precision blockin (6,ndima)    ! [input] matrix 
      double precision blockout(-2:2,ndima) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i
*rak:      integer ls     ! [fixed at 2] angular momentum of block
      double precision dt2m2, dt5m1, dt10, dt40, dt60, dt31, dt12, dt42
      double precision bin1, bin2, bin3, bin4, bin5, bin6
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=2
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_sd_a: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c  spint(-2) = cint(2)*dtrans(2,-2,2)
c  spint(-1) = cint(5)*dtrans(5,-1,2)
c  spint( 0) = cint(1)*dtrans(1,0,2) + cint(4)*dtrans(4,0,2) + cint(6)*dtrans(6,0,2)
c  spint( 1) = cint(3)*dtrans(3,1,2)
c  spint( 2) = cint(1)*dtrans(1,2,2) + cint(4)*dtrans(4,2,2)
      dt2m2 = dtrans(2,-2,2)
      dt5m1 = dtrans(5,-1,2)
      dt10  = dtrans(1, 0,2)
      dt40  = dtrans(4, 0,2)
      dt60  = dtrans(6, 0,2)
      dt31  = dtrans(3, 1,2)
      dt12  = dtrans(1, 2,2)
      dt42  = dtrans(4, 2,2)
      do 00100 i=1,ndima
        bin1 = blockin(1,i)
        bin2 = blockin(2,i)
        bin3 = blockin(3,i)
        bin4 = blockin(4,i)
        bin5 = blockin(5,i)
        bin6 = blockin(6,i)
        blockout(-2,i) = bin2*dt2m2
        blockout(-1,i) = bin5*dt5m1
        blockout( 0,i) = bin1*dt10 +
     &                   bin4*dt40 +
     &                   bin6*dt60
        blockout( 1,i) = bin3*dt31
        blockout( 2,i) = bin1*dt12 +
     &                   bin4*dt42
00100 continue
      end
*.......................................................................
      subroutine spcart_a_sd_b(blockin, blockout,
     &    ndima, ndimb, in_place)
      implicit none
c
c  transforms a block of integrals where the d function is the middle dimension
c  e.g., blockin(ndima,6,ndimb) and blockout(ndima,5,ndimb)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] trailing dimension of block
      double precision blockin (ndima, 6   ,ndimb) ! [input] matrix 
      double precision blockout(ndima, -2:2,ndimb) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i,j
*rak:      integer ls     ! [fixed at 2] angular momentum of block
      double precision dt2m2, dt5m1, dt10, dt40, dt60, dt31, dt12, dt42
      double precision bin1, bin2, bin3, bin4, bin5, bin6
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=2
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_sd_b: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c  spint(-2) = cint(2)*dtrans(2,-2,2)
c  spint(-1) = cint(5)*dtrans(5,-1,2)
c  spint( 0) = cint(1)*dtrans(1,0,2) + cint(4)*dtrans(4,0,2) + cint(6)*dtrans(6,0,2)
c  spint( 1) = cint(3)*dtrans(3,1,2)
c  spint( 2) = cint(1)*dtrans(1,2,2) + cint(4)*dtrans(4,2,2)
      dt2m2 = dtrans(2,-2,2)
      dt5m1 = dtrans(5,-1,2)
      dt10  = dtrans(1, 0,2)
      dt40  = dtrans(4, 0,2)
      dt60  = dtrans(6, 0,2)
      dt31  = dtrans(3, 1,2)
      dt12  = dtrans(1, 2,2)
      dt42  = dtrans(4, 2,2)
      do 00100 j=1,ndimb
        do 00200 i=1,ndima
          bin1 = blockin(i,1,j)
          bin2 = blockin(i,2,j)
          bin3 = blockin(i,3,j)
          bin4 = blockin(i,4,j)
          bin5 = blockin(i,5,j)
          bin6 = blockin(i,6,j)
          blockout(i,-2,j) = bin2*dt2m2
          blockout(i,-1,j) = bin5*dt5m1
          blockout(i, 0,j) = bin1*dt10 +
     &                       bin4*dt40 +
     &                       bin6*dt60
          blockout(i, 1,j) = bin3*dt31
          blockout(i, 2,j) = bin1*dt12 +
     &                       bin4*dt42
00200   continue
00100 continue
      end
*.......................................................................
      subroutine spcart_a_sf(blockin, blockout, ndima, in_place)
      implicit none
c
c  transforms a block of integrals where the f function is the slowest dimension
c  e.g., blockin(ndima,10) and blockout(ndima,7)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      double precision blockin (ndima,10)    ! [input] matrix 
      double precision blockout(ndima,-3:3)  ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i
*rak:      integer ls     ! [fixed at 3] angular momentum of block
      double precision dt2_m3, dt7_m3, dt5_m2, dt2_m1, dt7_m1, dt9_m1
      double precision dt3_0, dt8_0, dt10_0, dt1_1, dt4_1, dt6_1
      double precision dt3_2, dt8_2, dt1_3, dt4_3
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=3
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_sf:spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
*   spint(-3) = cint(2)*dtrans(2,-3,3) + cint(7)*dtrans(7,-3,3)
*   spint(-2) = cint(5)*dtrans(5,-2,3) 
*   spint(-1) = cint(2)*dtrans(2,-1,3) + cint(7)*dtrans(7,-1,3) + cint(9)*dtrans(9,-1,3)
*   spint( 0) = cint(3)*dtrans(3,0,3) + cint(8)*dtrans(8,0,3) + cint(10)*dtrans(10,0,3)
*   spint( 1) = cint(1)*dtrans(1,1,3) + cint(4)*dtrans(4,1,3) + cint(6)*dtrans(6,1,3)
*   spint( 2) = cint(3)*dtrans(3,2,3) + cint(8)*dtrans(8,2,3)
*   spint( 3) = cint(1)*dtrans(1,3,3) + cint(4)*dtrans(4,3,3)
      dt2_m3 = dtrans(2,-3,3)
      dt7_m3 = dtrans(7,-3,3)
      dt5_m2 = dtrans(5,-2,3) 
      dt2_m1 = dtrans(2,-1,3)
      dt7_m1 = dtrans(7,-1,3)
      dt9_m1 = dtrans(9,-1,3)
      dt3_0  = dtrans(3,0,3)
      dt8_0  = dtrans(8,0,3)
      dt10_0 = dtrans(10,0,3)
      dt1_1  = dtrans(1,1,3)
      dt4_1  = dtrans(4,1,3) 
      dt6_1  = dtrans(6,1,3)
      dt3_2  = dtrans(3,2,3)
      dt8_2  = dtrans(8,2,3)
      dt1_3  = dtrans(1,3,3)
      dt4_3  = dtrans(4,3,3)
c
      do 00100 i=1,ndima
        bin1 = blockin(i,1)
        bin2 = blockin(i,2)
        bin3 = blockin(i,3)
        bin4 = blockin(i,4)
        bin5 = blockin(i,5)
        bin6 = blockin(i,6)
        bin7 = blockin(i,7)
        bin8 = blockin(i,8)
        bin9 = blockin(i,9)
        bin10 = blockin(i,10)
        blockout(i,-3) = bin2*dt2_m3 + bin7*dt7_m3
        blockout(i,-2) = bin5*dt5_m2
        blockout(i,-1) = bin2*dt2_m1 + bin7*dt7_m1 + bin9*dt9_m1
        blockout(i, 0) = bin3*dt3_0 + bin8*dt8_0 + bin10*dt10_0
        blockout(i, 1) = bin1*dt1_1 + bin4*dt4_1 + bin6*dt6_1
        blockout(i, 2) = bin3*dt3_2 + bin8*dt8_2
        blockout(i, 3) = bin1*dt1_3 + bin4*dt4_3
00100 continue
      end
*.......................................................................
      subroutine spcart_sf_a(blockin, blockout, ndima, in_place)
      implicit none
c
c  transforms a block of integrals where the f function is the fastest dimension
c  e.g., blockin(10,ndima) and blockout(7,ndima)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      double precision blockin (10,ndima)    ! [input] matrix 
      double precision blockout(-3:3,ndima)  ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i
*rak:      integer ls     ! [fixed at 3] angular momentum of block
      double precision dt2_m3, dt7_m3, dt5_m2, dt2_m1, dt7_m1, dt9_m1
      double precision dt3_0, dt8_0, dt10_0, dt1_1, dt4_1, dt6_1
      double precision dt3_2, dt8_2, dt1_3, dt4_3
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=3
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_sf_a: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
*   spint(-3) = cint(2)*dtrans(2,-3,3) + cint(7)*dtrans(7,-3,3)
*   spint(-2) = cint(5)*dtrans(5,-2,3) 
*   spint(-1) = cint(2)*dtrans(2,-1,3) + cint(7)*dtrans(7,-1,3) + cint(9)*dtrans(9,-1,3)
*   spint( 0) = cint(3)*dtrans(3,0,3) + cint(8)*dtrans(8,0,3) + cint(10)*dtrans(10,0,3)
*   spint( 1) = cint(1)*dtrans(1,1,3) + cint(4)*dtrans(4,1,3) + cint(6)*dtrans(6,1,3)
*   spint( 2) = cint(3)*dtrans(3,2,3) + cint(8)*dtrans(8,2,3)
*   spint( 3) = cint(1)*dtrans(1,3,3) + cint(4)*dtrans(4,3,3)
      dt2_m3 = dtrans(2,-3,3)
      dt7_m3 = dtrans(7,-3,3)
      dt5_m2 = dtrans(5,-2,3) 
      dt2_m1 = dtrans(2,-1,3)
      dt7_m1 = dtrans(7,-1,3)
      dt9_m1 = dtrans(9,-1,3)
      dt3_0  = dtrans(3,0,3)
      dt8_0  = dtrans(8,0,3)
      dt10_0 = dtrans(10,0,3)
      dt1_1  = dtrans(1,1,3)
      dt4_1  = dtrans(4,1,3) 
      dt6_1  = dtrans(6,1,3)
      dt3_2  = dtrans(3,2,3)
      dt8_2  = dtrans(8,2,3)
      dt1_3  = dtrans(1,3,3)
      dt4_3  = dtrans(4,3,3)
c
      do 00100 i=1,ndima
        bin1  = blockin(1,i)
        bin2  = blockin(2,i)
        bin3  = blockin(3,i)
        bin4  = blockin(4,i)
        bin5  = blockin(5,i)
        bin6  = blockin(6,i)
        bin7  = blockin(7,i)
        bin8  = blockin(8,i)
        bin9  = blockin(9,i)
        bin10 = blockin(10,i)
        blockout(-3,i) = bin2*dt2_m3 + bin7*dt7_m3
        blockout(-2,i) = bin5*dt5_m2
        blockout(-1,i) = bin2*dt2_m1 + bin7*dt7_m1 + bin9*dt9_m1
        blockout( 0,i) = bin3*dt3_0 + bin8*dt8_0 + bin10*dt10_0
        blockout( 1,i) = bin1*dt1_1 + bin4*dt4_1 + bin6*dt6_1
        blockout( 2,i) = bin3*dt3_2 + bin8*dt8_2
        blockout( 3,i) = bin1*dt1_3 + bin4*dt4_3
00100 continue
      end
*.......................................................................
      subroutine spcart_a_sf_b(blockin, blockout,
     &    ndima, ndimb, in_place)
      implicit none
c
c  transforms a block of integrals where the f function is the middle dimension
c  e.g., blockin(ndima,10,ndimb) and blockout(ndima,7,ndimb)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] trailing dimension of block
      double precision blockin (ndima,10,ndimb)    ! [input] matrix 
      double precision blockout(ndima,-3:3,ndimb)  ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i,j
*rak:      integer ls     ! [fixed at 3] angular momentum of block
      double precision dt2_m3, dt7_m3, dt5_m2, dt2_m1, dt7_m1, dt9_m1
      double precision dt3_0, dt8_0, dt10_0, dt1_1, dt4_1, dt6_1
      double precision dt3_2, dt8_2, dt1_3, dt4_3
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=3
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_sf_b: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
*   spint(-3) = cint(2)*dtrans(2,-3,3) + cint(7)*dtrans(7,-3,3)
*   spint(-2) = cint(5)*dtrans(5,-2,3) 
*   spint(-1) = cint(2)*dtrans(2,-1,3) + cint(7)*dtrans(7,-1,3) + cint(9)*dtrans(9,-1,3)
*   spint( 0) = cint(3)*dtrans(3,0,3) + cint(8)*dtrans(8,0,3) + cint(10)*dtrans(10,0,3)
*   spint( 1) = cint(1)*dtrans(1,1,3) + cint(4)*dtrans(4,1,3) + cint(6)*dtrans(6,1,3)
*   spint( 2) = cint(3)*dtrans(3,2,3) + cint(8)*dtrans(8,2,3)
*   spint( 3) = cint(1)*dtrans(1,3,3) + cint(4)*dtrans(4,3,3)
      dt2_m3 = dtrans(2,-3,3)
      dt7_m3 = dtrans(7,-3,3)
      dt5_m2 = dtrans(5,-2,3) 
      dt2_m1 = dtrans(2,-1,3)
      dt7_m1 = dtrans(7,-1,3)
      dt9_m1 = dtrans(9,-1,3)
      dt3_0  = dtrans(3,0,3)
      dt8_0  = dtrans(8,0,3)
      dt10_0 = dtrans(10,0,3)
      dt1_1  = dtrans(1,1,3)
      dt4_1  = dtrans(4,1,3) 
      dt6_1  = dtrans(6,1,3)
      dt3_2  = dtrans(3,2,3)
      dt8_2  = dtrans(8,2,3)
      dt1_3  = dtrans(1,3,3)
      dt4_3  = dtrans(4,3,3)
c
      do 00100 j=1,ndimb
        do 00200 i=1,ndima
          bin1  = blockin(i,1,j)
          bin2  = blockin(i,2,j)
          bin3  = blockin(i,3,j)
          bin4  = blockin(i,4,j)
          bin5  = blockin(i,5,j)
          bin6  = blockin(i,6,j)
          bin7  = blockin(i,7,j)
          bin8  = blockin(i,8,j)
          bin9  = blockin(i,9,j)
          bin10 = blockin(i,10,j)
          blockout(i,-3,j) = bin2*dt2_m3 + bin7*dt7_m3
          blockout(i,-2,j) = bin5*dt5_m2
          blockout(i,-1,j) = bin2*dt2_m1 + bin7*dt7_m1 + bin9*dt9_m1
          blockout(i, 0,j) = bin3*dt3_0 + bin8*dt8_0 + bin10*dt10_0
          blockout(i, 1,j) = bin1*dt1_1 + bin4*dt4_1 + bin6*dt6_1
          blockout(i, 2,j) = bin3*dt3_2 + bin8*dt8_2
          blockout(i, 3,j) = bin1*dt1_3 + bin4*dt4_3
00200   continue
00100 continue
      end
*.......................................................................
      subroutine spcart_a_sg(blockin, blockout, ndima, in_place)
      implicit none
c
c  transforms a block of integrals where the g function is the slowest dimension
c  e.g., blockin(ndima,15) and blockout(ndima,9)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      double precision blockin (ndima,15)   ! [input] matrix 
      double precision blockout(ndima,-4:4) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i
      double precision dt2_m4, dt7_m4, dt5_m3, dt12_m3
      double precision dt2_m2, dt7_m2, dt9_m2, dt5_m1, dt12_m1, dt14_m1
      double precision dt1_0, dt4_0, dt6_0, dt11_0, dt13_0, dt15_0
      double precision dt3_1, dt8_1, dt10_1
      double precision dt1_2, dt6_2, dt11_2, dt13_2
      double precision dt3_3, dt8_3, dt1_4, dt4_4, dt11_4
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
      double precision bin11, bin12, bin13, bin14, bin15
*rak:      integer ls     ! [fixed at 4] angular momentum of block
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=4
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_sg: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
c  spint(-4) = cint(2)*dtrans(2,-4,4) + cint(7)*dtrans(7,-4,4)
c  spint(-3) = cint(5)*dtrans(5,-3,4) + cint(12)*dtrans(12,-3,4)
c  spint(-2) = cint(2)*dtrans(2,-2,4) + cint(7)*dtrans(7,-2,4) + cint(9)*dtrans(9,-2,4)
c  spint(-1) = cint(5)*dtrans(5,-1,4) + cint(12)*dtrans(12,-1,4) + cint(14)*dtrans(14,-1,4)
c  spint( 0) = cint(1)*dtrans(1,0,4)   + cint(4)*dtrans(4,0,4)   + cint(6)*dtrans(6,0,4)
c            + cint(11)*dtrans(11,0,4) + cint(13)*dtrans(13,0,4) + cint(15)*dtrans(15,0,4)
c  spint( 1) = cint(3)*dtrans(3,1,4) + cint(8)*dtrans(8,1,4) + cint(10)*dtrans(10,1,4)
c  spint( 2) = cint(1)*dtrans(1,2,4) + cint(6)*dtrans(6,2,4) + cint(11)*dtrans(11,2,4)
c            + cint(13)*dtrans(13,2,4)
c  spint( 3) = cint(3)*dtrans(3,3,4) + cint(8)*dtrans(8,3,4)
c  spint( 4) = cint(1)*dtrans(1,4,4) + cint(4)*dtrans(4,4,4) + cint(11)*dtrans(11,4,4)
c
      dt2_m4  = dtrans(2,-4,4)     
      dt7_m4  = dtrans(7,-4,4)     
      dt5_m3  = dtrans(5,-3,4)     
      dt12_m3 = dtrans(12,-3,4)    
      dt2_m2  = dtrans(2,-2,4)     
      dt7_m2  = dtrans(7,-2,4)     
      dt9_m2  = dtrans(9,-2,4)     
      dt5_m1  = dtrans(5,-1,4)     
      dt12_m1 = dtrans(12,-1,4)    
      dt14_m1 = dtrans(14,-1,4)    
      dt1_0   = dtrans(1,0,4)      
      dt4_0   = dtrans(4,0,4)      
      dt6_0   = dtrans(6,0,4)      
      dt11_0  = dtrans(11,0,4)     
      dt13_0  = dtrans(13,0,4)     
      dt15_0  = dtrans(15,0,4)     
      dt3_1   = dtrans(3,1,4)      
      dt8_1   = dtrans(8,1,4)      
      dt10_1  = dtrans(10,1,4)     
      dt1_2   = dtrans(1,2,4)      
      dt6_2   = dtrans(6,2,4)      
      dt11_2  = dtrans(11,2,4)     
      dt13_2  = dtrans(13,2,4)                 
      dt3_3   = dtrans(3,3,4)      
      dt8_3   = dtrans(8,3,4)      
      dt1_4   = dtrans(1,4,4)      
      dt4_4   = dtrans(4,4,4)      
      dt11_4  = dtrans(11,4,4)     
      do 00100 i=1,ndima
        bin1 = blockin(i,1)
        bin2 = blockin(i,2)
        bin3 = blockin(i,3)
        bin4 = blockin(i,4)
        bin5 = blockin(i,5)
        bin6 = blockin(i,6)
        bin7 = blockin(i,7)
        bin8 = blockin(i,8)
        bin9 = blockin(i,9)
        bin10 = blockin(i,10)
        bin11 = blockin(i,11)
        bin12 = blockin(i,12)
        bin13 = blockin(i,13)
        bin14 = blockin(i,14)
        bin15 = blockin(i,15)
c
        blockout(i,-4) = bin2*dt2_m4  + bin7*dt7_m4  
        blockout(i,-3) = bin5*dt5_m3  + bin12*dt12_m3 
        blockout(i,-2) = bin2*dt2_m2  + bin7*dt7_m2   + bin9*dt9_m2
        blockout(i,-1) = bin5*dt5_m1  + bin12*dt12_m1 + bin14*dt14_m1
        blockout(i, 0) = bin1*dt1_0   + bin4*dt4_0    + bin6*dt6_0   +
     &                   bin11*dt11_0 + bin13*dt13_0  + bin15*dt15_0
        blockout(i, 1) = bin3*dt3_1   + bin8*dt8_1    + bin10*dt10_1
        blockout(i, 2) = bin1*dt1_2   + bin6*dt6_2    + bin11*dt11_2 +
     &                   bin13*dt13_2
        blockout(i, 3) = bin3*dt3_3   + bin8*dt8_3
        blockout(i, 4) = bin1*dt1_4   + bin4*dt4_4    + bin11*dt11_4
00100 continue
      end
*.......................................................................
      subroutine spcart_sg_a(blockin, blockout, ndima, in_place)
      implicit none
c
c  transforms a block of integrals where the g function is the fastest dimension
c  e.g., blockin(15,ndima) and blockout(9,ndima)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      double precision blockin (15,ndima)   ! [input] matrix 
      double precision blockout(-4:4,ndima) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i
      double precision dt2_m4, dt7_m4, dt5_m3, dt12_m3
      double precision dt2_m2, dt7_m2, dt9_m2, dt5_m1, dt12_m1, dt14_m1
      double precision dt1_0, dt4_0, dt6_0, dt11_0, dt13_0, dt15_0
      double precision dt3_1, dt8_1, dt10_1
      double precision dt1_2, dt6_2, dt11_2, dt13_2
      double precision dt3_3, dt8_3, dt1_4, dt4_4, dt11_4
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
      double precision bin11, bin12, bin13, bin14, bin15
*rak:      integer ls     ! [fixed at 4] angular momentum of block
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=4
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_sg_a: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
c  spint(-4) = cint(2)*dtrans(2,-4,4) + cint(7)*dtrans(7,-4,4)
c  spint(-3) = cint(5)*dtrans(5,-3,4) + cint(12)*dtrans(12,-3,4)
c  spint(-2) = cint(2)*dtrans(2,-2,4) + cint(7)*dtrans(7,-2,4) + cint(9)*dtrans(9,-2,4)
c  spint(-1) = cint(5)*dtrans(5,-1,4) + cint(12)*dtrans(12,-1,4) + cint(14)*dtrans(14,-1,4)
c  spint( 0) = cint(1)*dtrans(1,0,4)   + cint(4)*dtrans(4,0,4)   + cint(6)*dtrans(6,0,4)
c            + cint(11)*dtrans(11,0,4) + cint(13)*dtrans(13,0,4) + cint(15)*dtrans(15,0,4)
c  spint( 1) = cint(3)*dtrans(3,1,4) + cint(8)*dtrans(8,1,4) + cint(10)*dtrans(10,1,4)
c  spint( 2) = cint(1)*dtrans(1,2,4) + cint(6)*dtrans(6,2,4) + cint(11)*dtrans(11,2,4)
c            + cint(13)*dtrans(13,2,4)
c  spint( 3) = cint(3)*dtrans(3,3,4) + cint(8)*dtrans(8,3,4)
c  spint( 4) = cint(1)*dtrans(1,4,4) + cint(4)*dtrans(4,4,4) + cint(11)*dtrans(11,4,4)
c
      dt2_m4  = dtrans(2,-4,4)     
      dt7_m4  = dtrans(7,-4,4)     
      dt5_m3  = dtrans(5,-3,4)     
      dt12_m3 = dtrans(12,-3,4)    
      dt2_m2  = dtrans(2,-2,4)     
      dt7_m2  = dtrans(7,-2,4)     
      dt9_m2  = dtrans(9,-2,4)     
      dt5_m1  = dtrans(5,-1,4)     
      dt12_m1 = dtrans(12,-1,4)    
      dt14_m1 = dtrans(14,-1,4)    
      dt1_0   = dtrans(1,0,4)      
      dt4_0   = dtrans(4,0,4)      
      dt6_0   = dtrans(6,0,4)      
      dt11_0  = dtrans(11,0,4)     
      dt13_0  = dtrans(13,0,4)     
      dt15_0  = dtrans(15,0,4)     
      dt3_1   = dtrans(3,1,4)      
      dt8_1   = dtrans(8,1,4)      
      dt10_1  = dtrans(10,1,4)     
      dt1_2   = dtrans(1,2,4)      
      dt6_2   = dtrans(6,2,4)      
      dt11_2  = dtrans(11,2,4)     
      dt13_2  = dtrans(13,2,4)                 
      dt3_3   = dtrans(3,3,4)      
      dt8_3   = dtrans(8,3,4)      
      dt1_4   = dtrans(1,4,4)      
      dt4_4   = dtrans(4,4,4)      
      dt11_4  = dtrans(11,4,4)     
      do 00100 i=1,ndima
        bin1  = blockin( 1,i)
        bin2  = blockin( 2,i)
        bin3  = blockin( 3,i)
        bin4  = blockin( 4,i)
        bin5  = blockin( 5,i)
        bin6  = blockin( 6,i)
        bin7  = blockin( 7,i)
        bin8  = blockin( 8,i)
        bin9  = blockin( 9,i)
        bin10 = blockin(10,i)
        bin11 = blockin(11,i)
        bin12 = blockin(12,i)
        bin13 = blockin(13,i)
        bin14 = blockin(14,i)
        bin15 = blockin(15,i)
c
        blockout(-4,i) = bin2*dt2_m4  + bin7*dt7_m4  
        blockout(-3,i) = bin5*dt5_m3  + bin12*dt12_m3 
        blockout(-2,i) = bin2*dt2_m2  + bin7*dt7_m2   + bin9*dt9_m2
        blockout(-1,i) = bin5*dt5_m1  + bin12*dt12_m1 + bin14*dt14_m1
        blockout( 0,i) = bin1*dt1_0   + bin4*dt4_0    + bin6*dt6_0   +
     &                   bin11*dt11_0 + bin13*dt13_0  + bin15*dt15_0
        blockout( 1,i) = bin3*dt3_1   + bin8*dt8_1    + bin10*dt10_1
        blockout( 2,i) = bin1*dt1_2   + bin6*dt6_2    + bin11*dt11_2 +
     &                   bin13*dt13_2
        blockout( 3,i) = bin3*dt3_3   + bin8*dt8_3
        blockout( 4,i) = bin1*dt1_4   + bin4*dt4_4    + bin11*dt11_4
00100 continue
      end
*.......................................................................
      subroutine spcart_a_sg_b(blockin, blockout,
     &    ndima, ndimb, in_place)
      implicit none
c
c  transforms a block of integrals where the g function is the middle dimension
c  e.g., blockin(ndima,15,ndimb) and blockout(ndima,9,ndimb)
c
#include "mafdecls.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] trailing dimension of block
      double precision blockin (ndima,15,ndimb)   ! [input] matrix 
      double precision blockout(ndima,-4:4,ndimb) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
c: local
      integer i,j
      double precision dt2_m4, dt7_m4, dt5_m3, dt12_m3
      double precision dt2_m2, dt7_m2, dt9_m2, dt5_m1, dt12_m1, dt14_m1
      double precision dt1_0, dt4_0, dt6_0, dt11_0, dt13_0, dt15_0
      double precision dt3_1, dt8_1, dt10_1
      double precision dt1_2, dt6_2, dt11_2, dt13_2
      double precision dt3_3, dt8_3, dt1_4, dt4_4, dt11_4
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
      double precision bin11, bin12, bin13, bin14, bin15
*rak:      integer ls     ! [fixed at 4] angular momentum of block
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=4
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_sg_b: spcart not initialized properly',
     &    sph_cart_init)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
c  spint(-4) = cint(2)*dtrans(2,-4,4) + cint(7)*dtrans(7,-4,4)
c  spint(-3) = cint(5)*dtrans(5,-3,4) + cint(12)*dtrans(12,-3,4)
c  spint(-2) = cint(2)*dtrans(2,-2,4) + cint(7)*dtrans(7,-2,4) + cint(9)*dtrans(9,-2,4)
c  spint(-1) = cint(5)*dtrans(5,-1,4) + cint(12)*dtrans(12,-1,4) + cint(14)*dtrans(14,-1,4)
c  spint( 0) = cint(1)*dtrans(1,0,4)   + cint(4)*dtrans(4,0,4)   + cint(6)*dtrans(6,0,4)
c            + cint(11)*dtrans(11,0,4) + cint(13)*dtrans(13,0,4) + cint(15)*dtrans(15,0,4)
c  spint( 1) = cint(3)*dtrans(3,1,4) + cint(8)*dtrans(8,1,4) + cint(10)*dtrans(10,1,4)
c  spint( 2) = cint(1)*dtrans(1,2,4) + cint(6)*dtrans(6,2,4) + cint(11)*dtrans(11,2,4)
c            + cint(13)*dtrans(13,2,4)
c  spint( 3) = cint(3)*dtrans(3,3,4) + cint(8)*dtrans(8,3,4)
c  spint( 4) = cint(1)*dtrans(1,4,4) + cint(4)*dtrans(4,4,4) + cint(11)*dtrans(11,4,4)
c
      dt2_m4  = dtrans(2,-4,4)     
      dt7_m4  = dtrans(7,-4,4)     
      dt5_m3  = dtrans(5,-3,4)     
      dt12_m3 = dtrans(12,-3,4)    
      dt2_m2  = dtrans(2,-2,4)     
      dt7_m2  = dtrans(7,-2,4)     
      dt9_m2  = dtrans(9,-2,4)     
      dt5_m1  = dtrans(5,-1,4)     
      dt12_m1 = dtrans(12,-1,4)    
      dt14_m1 = dtrans(14,-1,4)    
      dt1_0   = dtrans(1,0,4)      
      dt4_0   = dtrans(4,0,4)      
      dt6_0   = dtrans(6,0,4)      
      dt11_0  = dtrans(11,0,4)     
      dt13_0  = dtrans(13,0,4)     
      dt15_0  = dtrans(15,0,4)     
      dt3_1   = dtrans(3,1,4)      
      dt8_1   = dtrans(8,1,4)      
      dt10_1  = dtrans(10,1,4)     
      dt1_2   = dtrans(1,2,4)      
      dt6_2   = dtrans(6,2,4)      
      dt11_2  = dtrans(11,2,4)     
      dt13_2  = dtrans(13,2,4)                 
      dt3_3   = dtrans(3,3,4)      
      dt8_3   = dtrans(8,3,4)      
      dt1_4   = dtrans(1,4,4)      
      dt4_4   = dtrans(4,4,4)      
      dt11_4  = dtrans(11,4,4)     
      do 00100 j=1,ndimb
        do 00200 i=1,ndima
          bin1  = blockin(i, 1,j)
          bin2  = blockin(i, 2,j)
          bin3  = blockin(i, 3,j)
          bin4  = blockin(i, 4,j)
          bin5  = blockin(i, 5,j)
          bin6  = blockin(i, 6,j)
          bin7  = blockin(i, 7,j)
          bin8  = blockin(i, 8,j)
          bin9  = blockin(i, 9,j)
          bin10 = blockin(i,10,j)
          bin11 = blockin(i,11,j)
          bin12 = blockin(i,12,j)
          bin13 = blockin(i,13,j)
          bin14 = blockin(i,14,j)
          bin15 = blockin(i,15,j)
c
          blockout(i,-4,j) = bin2*dt2_m4  + bin7*dt7_m4  
          blockout(i,-3,j) = bin5*dt5_m3  + bin12*dt12_m3 
          blockout(i,-2,j) = bin2*dt2_m2  + bin7*dt7_m2   +
     &                       bin9*dt9_m2
          blockout(i,-1,j) = bin5*dt5_m1  + bin12*dt12_m1 +
     &                       bin14*dt14_m1
          blockout(i, 0,j) = bin1*dt1_0   + bin4*dt4_0    +
     &                       bin6*dt6_0   + bin11*dt11_0  +
     &                       bin13*dt13_0 + bin15*dt15_0
          blockout(i, 1,j) = bin3*dt3_1   + bin8*dt8_1    +
     &                       bin10*dt10_1
          blockout(i, 2,j) = bin1*dt1_2   + bin6*dt6_2    +
     &                       bin11*dt11_2 + bin13*dt13_2
          blockout(i, 3,j) = bin3*dt3_3   + bin8*dt8_3
          blockout(i, 4,j) = bin1*dt1_4   + bin4*dt4_4    +
     &                       bin11*dt11_4
00200   continue
00100 continue
      end
*.......................................................................
      subroutine spcart_tran1e(
     &    buf, scr,
     &    nbf_xr,nbf_xc,type_r,
     &    nbf_sr,nbf_sc,type_c,
     &    print)
      implicit none
#include "stdio.fh"
c
c routine that transforms a 1e cartesian block buf_cart(nbf_xr,nbf_xc) to 
c    a spherical block buf_sph(nbf_sr,nbf_sc) 
c  
c  x --> implies cartesian
c  s --> implies spherical
c  r --> implies row
c  c --> implies column
c 
c  remember that a Ov block for ish and jsh is 
c      from the integral api Ov(jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c      
      integer nbf_xr, nbf_xc  ! [input] size of cartesian block
      integer nbf_sr, nbf_sc  ! [input] size of spherical block
      integer type_r, type_c  ! [input] angular momentem for r and c
      double precision buf(*) ! [input/output] cartesian block on input
*.............................!   and spherical block on output      
      double precision scr(*) ! [scratch] use to hold half transformed block
      logical print           ! [input} print integrals at each stage of the
*.............................!   transformation (cart/half/spherical)
c
      if (type_r.lt.2.and.type_c.lt.2) then
c.................................. neither c or r need to be transformed 
c                                   (X,Y) X is s, p, or l and Y is s, p, or l
        if (print) then
          write(luout,*)
     &        ' cartesian matrix and spherical matrix the same '
          call output(buf,1,nbf_xr,1,nbf_xc,nbf_xr,nbf_xc,1)
        endif
c
      elseif (type_r.lt.2.and.type_c.ge.2) then
c.................................. c needs to be transformed
* print cartesian matrix  
* buf is buf(spherical,cartesian)
        if (print) then
          write(luout,*)' cartesian matrix '
          call output(buf,1,nbf_xr,1,nbf_xc,nbf_xr,nbf_xc,1)
        endif
* scr is buf(spherical,cartesian) ! copy it
        call dcopy((nbf_xr*nbf_xc),buf,1,scr,1)
        if (nbf_xr.ne.nbf_sr) call errquit
     &      ('spcart_tran1e: nbf_xr.ne.nbf_sr  (xr-sr) =',
     &      (nbf_xr-nbf_sr))
        call spcart_a_s(scr,buf,nbf_sr,type_c,.false.)
        if (print) then
          write(luout,*)' spherical matrix '
          call output(buf,1,nbf_sr,1,nbf_sc,nbf_sr,nbf_sc,1)
        endif
      elseif (type_r.ge.2.and.type_c.lt.2) then
c.................................. r needs to be transformed
* print cartesian matrix  
* buf is buf(cartesian,cartesian)
        if (print) then
          write(luout,*)' cartesian matrix '
          call output(buf,1,nbf_xr,1,nbf_xc,nbf_xr,nbf_xc,1)
        endif
* scr is buf(spherical,cartesian) ! copy it
        call dcopy((nbf_xr*nbf_xc),buf,1,scr,1)
        if (nbf_xc.ne.nbf_sc) call errquit
     &      ('spcart_tran1e: nbf_xc.ne.nbf_sc  (xc-sc) =',
     &      (nbf_xc-nbf_sc))
        call spcart_s_a(scr,buf,nbf_sc,type_r,.false.)
        if (print) then
          write(luout,*)' spherical matrix '
          call output(buf,1,nbf_sr,1,nbf_sc,nbf_sr,nbf_sc,1)
        endif
      elseif (type_r.ge.2.and.type_c.ge.2) then
c.................................. both r and c need to be transformed
* print cartesian matrix  
* buf is buf(cartesian,cartesian)
        if (print) then
          write(luout,*)' cartesian matrix '
          call output(buf,1,nbf_xr,1,nbf_xc,nbf_xr,nbf_xc,1)
        endif
*
*... buf(xr,xc) -> scr(xr,sc) : scr is half transformed matrix
        call spcart_a_s(buf,scr,nbf_xr,type_c,.false.)
* print half transformed matrix
        if (print) then
          write(luout,*)' half cartesian half spherical matrix '
          call output(scr,1,nbf_xr,1,nbf_sc,nbf_xr,nbf_sc,1)
        endif
*
*... scr(xr,sc) -> buf(sr,sc)
        call spcart_s_a(scr,buf,nbf_sc,type_r,.false.)
* print spherical block
        if (print) then
          write(luout,*)' half cartesian half spherical matrix '
          call output(buf,1,nbf_sr,1,nbf_sc,nbf_sr,nbf_sc,1)
        endif
      else
        write(luout,*) ' case not possible '
        call errquit('spcart_tran1e: should never get here ?? ',911)
      endif
c
      end
*.......................................................................
      subroutine spcart_bra2etran(
     &    buf, scr,
     &    nbf_xj,nbf_xi,
     &    nbf_sj,nbf_si,
     &    type_j,type_i,
     &    ndim_ket,
     &    print)
      implicit none
#include "stdio.fh"
c
c routine that transforms a 2e cartesian block buf_cart(ndim_ket,nbf_xr,nbf_xc) to 
c    a spherical block buf_sph(ndim_ket,nbf_sr,nbf_sc) 
c  
c  x --> implies cartesian
c  s --> implies spherical
c 
c  remember that a 2e block for ish jsh for any k/l sh is 
c      from the integral api ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c
c row-j col-i
c
c  blockin(ndim_ket,jxR,ixR)*trans = blockout(ndim_ket,jsR,isR)
c      
c
      integer nbf_xj, nbf_xi  ! [input] size of cartesian block
      integer nbf_sj, nbf_si  ! [input] size of spherical block
      integer type_j, type_i  ! [input] angular momentem for j and i
      integer ndim_ket
      double precision buf(*) ! [input/output] cartesian block on input
*.............................!   and spherical block on output      
      double precision scr(*) ! [scratch] use to hold half transformed block
      logical print           ! [input} print integrals at each stage of the
*.............................!   transformation (cart/half/spherical)
*::local
      integer count, iket, jr, ir
c
      if (type_j.lt.2.and.type_i.lt.2) then
c.................................. neither i or j need to be transformed 
c                                   (X,Y) X is s, p, or l and Y is s, p, or l
        if (print) then
          write(luout,*)
     &        ' cartesian matrix and spherical matrix the same '
          count = 1
          do ir = 1,nbf_xi
            do jr = 1,nbf_xj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
*.......make this pretty later
*rak:c  buf has integrals
*rak:c  scr has slices of buf
*rak:        if (print) then
*rak:          write(luout,*)
*rak:     &        ' cartesian matrix and spherical matrix the same '
*rak:          do iket = 0,(ndim_ket-1)
*rak:            scr_count = 0
*rak:            do cc = 1,nbf_xi 
*rak:              do rr = 1,nbf_xj
*rak:                offset = ndim_ket*scr_count + 1
*rak:                scr_count = scr_count + 1
*rak:                scr(scr_count) = buf(offset)
*rak:              enddo
*rak:            enddo
*rak:          enddo
*rak:          call output(scr,1,nbf_xj,1,nbf_xi,nbf_xj,nbf_xi,1)
*rak:        endif
      elseif (type_j.lt.2.and.type_i.ge.2) then
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c.................................. i needs to be transformed
* print cartesian matrix  
* buf is buf(ndim_ket,j_spherical,i_cartesian)
        if (print) then
          write(luout,*)' (ket:spherical:cartesian) matrix '
          count = 1
          do ir = 1,nbf_xi
            do jr = 1,nbf_sj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
* scr is buf(ndim_ket,j_spherical,i_cartesian) ! copy it
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
        call dcopy((ndim_ket*nbf_sj*nbf_xi),buf,1,scr,1)
        if (nbf_xj.ne.nbf_sj) call errquit
     &      ('spcart_bra2etran: nbf_xj.ne.nbf_sj  (xj-sj) =',
     &      (nbf_xj-nbf_sj))
        call spcart_a_s(scr,buf,(ndim_ket*nbf_sj),type_i,.false.)
        if (print) then
          write(luout,*)' (ket:spherical:shperical) matrix '
          count = 1
          do ir = 1,nbf_si
            do jr = 1,nbf_sj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
      elseif (type_j.ge.2.and.type_i.lt.2) then
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c.................................. j needs to be transformed
* print cartesian matrix  
* buf is buf(ndim_ket,j_cartesian,i_spherical)
        if (print) then
          write(luout,*)' (ket:cartesian:spherical) matrix '
          count = 1
          do ir = 1,nbf_si
            do jr = 1,nbf_xj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
* buf is buf(ndim_ket,j_cartesian,i_spherical)
* scr is buf(ndim_ket,j_cartesian,i_spherical) ! copy it
        call dcopy((ndim_ket*nbf_xj*nbf_si),buf,1,scr,1)
        if (nbf_xi.ne.nbf_si) call errquit
     &      ('spcart_bra2etran: nbf_xc.ne.nbf_sc  (xi-si) =',
     &      (nbf_xi-nbf_si))
        call spcart_a_s_b(scr,buf,ndim_ket,nbf_si,type_j,.false.)
        if (print) then
          write(luout,*)' (ket:spherical:spherical) matrix '
          count = 1
          do ir = 1,nbf_si
            do jr = 1,nbf_sj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
      elseif (type_j.ge.2.and.type_i.ge.2) then
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c........................ both j and i need to be transformed
* print cartesian matrix  
* buf is buf(ndim_ket,j_cartesian,i_cartesian)
        if (print) then
          write(luout,*)' (ket:cartesian:cartesian) matrix '
          count = 1
          do ir = 1,nbf_xi
            do jr = 1,nbf_xj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
*
*... buf(ndim_ket,xj,xi) -> scr(xj,si) : scr is half transformed matrix
        call spcart_a_s(buf,scr,(ndim_ket*nbf_xj),type_i,.false.)
* print half transformed matrix
        if (print) then
          write(luout,*)' (ket:cartesian:spherical) matrix '
          count = 1
          do ir = 1,nbf_si
            do jr = 1,nbf_xj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,scr(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
*
*... scr(xj,si) -> buf(sj,si)
        call spcart_a_s_b(scr,buf,ndim_ket,nbf_si,type_j,.false.)
* print spherical block
        if (print) then
          write(luout,*)' (ket:spherical:spherical) matrix '
          count = 1
          do ir = 1,nbf_si
            do jr = 1,nbf_sj
              do iket = 1,ndim_ket
                write(luout,10000)iket,jr,ir,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
c
      else
        write(luout,*) ' case not possible '
        call errquit('spcart_bra2etran: should never get here ?? ',911)
      endif
c
10000 format('<ket=',i7,'|j=',i6,'|i=',i6,':count=',i10,
     &    '> =',1pd20.10)
      end
*.......................................................................
      subroutine spcart_ket2etran(
     &    buf, scr,
     &    nbf_xl,nbf_xk,
     &    nbf_sl,nbf_sk,
     &    type_l,type_k,
     &    ndim_bra,
     &    print)
      implicit none
#include "stdio.fh"
c
c routine that transforms a 2e cartesian block 
c     buf_cart(nbf_xl,nbf_xk,ndim_bra) to 
c    a spherical block buf_sph(nbf_sl,nbf_sk,ndim_bra) 
c  
c  x --> implies cartesian
c  s --> implies spherical
c 
c  remember that a 2e block for ksh lsh for any i(j)sh is 
c      from the integral api ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
c
c row-l col-k
c
c  blockin(lxR,kxR,ndim_bra)*trans = blockout(lsR,ksR,ndim_bra)
c      
c
      integer nbf_xl, nbf_xk  ! [input] size of cartesian block
      integer nbf_sl, nbf_sk  ! [input] size of spherical block
      integer type_l, type_k  ! [input] angular momentem for k and l
      integer ndim_bra
      double precision buf(*) ! [input/output] cartesian block on input
*.............................!   and spherical block on output      
      double precision scr(*) ! [scratch] use to hold half transformed block
      logical print           ! [input} print integrals at each stage of the
*.............................!   transformation (cart/half/spherical)
*::local
      integer count, ibra, lr, kr
c
      if (type_l.lt.2.and.type_k.lt.2) then
c.................................. neither k or l need to be transformed 
c                                   (X,Y) X is s, p, or l and Y is s, p, or l
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
        if (print) then
          write(luout,*)
     &        ' cartesian matrix and spherical matrix the same '
          count = 1
          do lr = 1,nbf_xl
            do kr = 1,nbf_xk
              do ibra = 1,ndim_bra
                write(luout,10000)ibra,lr,kr,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
      elseif (type_l.lt.2.and.type_k.ge.2) then
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
c.................................. k needs to be transformed
* print cartesian matrix  
* buf is buf(spherical,cartesian,ndim_bra)
        if (print) then
          write(luout,*)' (spherical:cartesian:bra) matrix '
          count = 1
          do ibra = 1,ndim_bra
            do kr = 1,nbf_xk
              do lr = 1,nbf_sl
                write(luout,10000)ibra,lr,kr,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
* scr is buf(l_spherical,k_cartesian,ndim_bra) ! copy it
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
        call dcopy((ndim_bra*nbf_sl*nbf_xk),buf,1,scr,1)
        if (nbf_xl.ne.nbf_sl) call errquit
     &      ('spcart_ket2etran: nbf_xl.ne.nbf_sl  (xl-sl) =',
     &      (nbf_xl-nbf_sl))
        call spcart_a_s_b(scr,buf,nbf_sl,ndim_bra,type_k,.false.)
        if (print) then
          write(luout,*)' (l-spherical:k-shperical:bra) matrix '
          count = 1
          do ibra = 1,ndim_bra
            do kr = 1,nbf_sk
              do lr = 1,nbf_sl
                write(luout,10000)ibra,lr,kr,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
      elseif (type_l.ge.2.and.type_k.lt.2) then
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
c.................................. l needs to be transformed
* print cartesian matrix  
* buf is buf(l_cartesian,k_spherical,ndim_bra)
        if (print) then
          write(luout,*)' (l-cartesian:k-spherical:bra) matrix '
          count = 1
          do ibra = 1,ndim_bra
            do lr = 1,nbf_xl
              do kr = 1,nbf_sk
                write(luout,10000)ibra,lr,kr,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
* buf is buf(l_cartesian,k_spherical,ndim_bra)
* scr is buf(l_cartesian,k_spherical,ndim_bra)
        call dcopy((ndim_bra*nbf_xl*nbf_sk),buf,1,scr,1)
        if (nbf_xk.ne.nbf_sk) call errquit
     &      ('spcart_ket2etran: nbf_xc.ne.nbf_sc  (xi-si) =',
     &      (nbf_xk-nbf_sk))
        call spcart_s_a(scr,buf,(ndim_bra*nbf_sk),type_l,.false.)
        if (print) then
          write(luout,*)' (l-spherical:k-spherical:bra) matrix '
          count = 1
          do ibra = 1,ndim_bra
            do lr = 1,nbf_sl
              do kr = 1,nbf_sk
                write(luout,10000)ibra,lr,kr,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
      elseif (type_l.ge.2.and.type_k.ge.2) then
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
c........................ both k and l need to be transformed
* print cartesian matrix  
* buf is buf(k_cartesian,l_cartesian,ndim_bra)
        if (print) then
          write(luout,*)' (l-cartesian:k-cartesian:bra) matrix '
          count = 1
          do ibra = 1,ndim_bra
            do lr = 1,nbf_xl
              do kr = 1,nbf_xk
                write(luout,10000)ibra,lr,kr,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
*
*... buf(xl,xk) -> scr(xl,sk,ndim_bra,) : scr is half transformed matrix
        call spcart_a_s_b(buf,scr,ndim_bra,nbf_xl,type_k,.false.)
* print half transformed matrix
        if (print) then
          write(luout,*)' (ket:cartesian:spherical) matrix '
          count = 1
          do ibra = 1,ndim_bra
            do lr = 1,nbf_xl
              do kr = 1,nbf_sk
                write(luout,10000)ibra,lr,kr,count,scr(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
*
*... scr(xl,sk) -> buf(sl,sk)
        call spcart_a_s_b(scr,buf,ndim_bra,nbf_sk,type_l,.false.)
* print spherical block
        if (print) then
          write(luout,*)' (ket:spherical:spherical) matrix '
          count = 1
          do ibra = 1,ndim_bra
            do lr = 1,nbf_sl
              do kr = 1,nbf_sk
                write(luout,10000)ibra,lr,kr,count,buf(count)
                count = count + 1
              enddo
            enddo
          enddo
        endif
c
      else
        write(luout,*) ' case not possible '
        call errquit('spcart_ket2etran: should never get here ?? ',911)
      endif
c
10000 format('<bra=',i7,'|l=',i6,'|k=',i6,':count=',i10,
     &    '> =',1pd20.10)
      end
