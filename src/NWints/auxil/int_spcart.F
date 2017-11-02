      logical function spcart_init(lmaxin,normalize,all_spherical)
c $Id$
*::cr::7
*--------------------------------------------------*
* COPYRIGHT (C) 1994, 1995, 1996, 1997, 1998, 1999 *
*         Pacific Northwest National Laboratory,   * 
*         Battelle Memorial Institute.             *
*--------------------------------------------------*
*------------> All Rights Reserved <---------------*
*--------------------------------------------------*
* 
* initialization of spherical to cartesian tranformation array
*... stored on heap.  
*... stored up to lmax values
*...
*   spcart(iccart,icsp,l)  => 
*                          spcart((lmax+1)*(lmax+2)/2,1:2*lmax+1,0:lmax)
*    lmax = 5 h functions  => size = 21*11*6 = 1386 
*   store array 34% sparse for simplicity.  1386 doubles is 11 Kbytes.
*
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "spcartP.fh"
c::passed
      integer lmaxin    ! [input] init transformed values up to lmaxin
      logical normalize ! [input] normalize the coefficients for 
*                                 integral transformations
*                         true for integral transformations
*                         false for ECP integral computations
*
      logical all_spherical ! [input] generate all spherical components 
*                                     e.g., 6 d sphericals etc.
*
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
     &      ('spcart_init: all spherical components not working yet',
     &      911, UNKNOWN_ERR)
c
      if (sph_cart_init.eq.SPH_CART_INIT_VALUE) then
        if (lmaxin.gt.lmax_init) then
          if (.not.spcart_terminate()) call errquit
     &          ('spcart_init: error terminating old spcart_init',911, 
     &          UNKNOWN_ERR)
        else
          spcart_init = .true.
          return           ! initialization already done to cover lmaxin
        endif
      endif
*
      if (all_spherical) then
        sph_cart_allsph = .true.
      else
        sph_cart_allsph = .false.
      endif
*
      call defNxyz(lmaxin)
c
      size_sp2c = lmaxin+1
      size_sp2c = size_sp2c*(2*lmaxin+1)
      size_sp2c = size_sp2c*(((lmaxin+1)*(lmaxin+2))/2)
c
      active_sp2c = 
     &      ma_alloc_get(mt_dbl,size_sp2c,'sph 2 cart trans array',
     &      h_sp2c,k_sp2c)
      if (.not.active_sp2c)  call errquit
     &      ('spcart_init: alloc_get failed for size',size_sp2c,
     &       MEM_ERR)
      call dcopy(size_sp2c,0.0d00,0,dbl_mb(k_sp2c),1)
c
* generate transformation matricies by recursion
      call xlmcoeff(lmaxin,dbl_mb(k_sp2c),normalize)
* generate all spherical components
      call xlm_coeff_full(lmaxin,dbl_mb(k_sp2c),normalize)
* allocate memory for index array
      active_sp2c_lindx = ma_alloc_get(
     &      mt_int,(lmaxin+1),' ptrs array xlm sph 2 cart ',
     &      h_sp2c_lindx,k_sp2c_lindx)
      if (.not.active_sp2c_lindx) call errquit
     &      ('spcart_init: alloc_get failed (index) ',911, MEM_ERR)
      call ifill((lmaxin+1),0,int_mb(k_sp2c_lindx),1)
* allocate memory for index array for inverse transform
      active_invsp2c_lindx = ma_alloc_get(
     &      mt_int,(lmaxin+1),
     &      ' ptrs array inverse xlm sph 2 cart ',
     &      h_invsp2c_lindx,k_invsp2c_lindx)
      if (.not.active_invsp2c_lindx) call errquit
     &      ('spcart_init: alloc_get failed (index) ',911, MEM_ERR)
      call ifill((lmaxin+1),0,int_mb(k_invsp2c_lindx),1)
* determine size of compressed transformation arrays
      l_block_size = 0
      do 00100 lval=0,Lmaxin
        l2 = (((lval+1)*(lval+2))/2)
        ls = (2*lval+1) 
        l_block_size = l_block_size + l2*ls
00100 continue
* allocate memory for compressed transformation arrays
      active_sp2c_cmp = ma_alloc_get
     &      (mt_dbl,l_block_size,'sph 2 cart trans array cmp',
     &      h_sp2c_cmp,k_sp2c_cmp)
      if (.not. active_sp2c_cmp) call errquit
     &      ('spcart_init: alloc_get failed (array cmp) ',911, MEM_ERR)
      call dcopy(l_block_size,0.0d00,0,dbl_mb(k_sp2c_cmp),1)
* allocate memory for compressed inverse transformation arrays
      active_invsp2c_cmp = ma_alloc_get
     &      (mt_dbl,l_block_size,'sph 2 cart trans array cmp',
     &      h_invsp2c_cmp,k_invsp2c_cmp)
      if (.not. active_invsp2c_cmp) call errquit
     &      ('spcart_init: alloc_get failed (inv array cmp) ',911,
     &      MEM_ERR)
      call dcopy(l_block_size,0.0d00,0,dbl_mb(k_invsp2c_cmp),1)
* allocate memory for scale vector
      active_cart_norm_scale= ma_alloc_get
     &      (mt_dbl,(2*lmaxin+1)*(lmaxin+1),
     &      'sph 2 cart scale arryy ',
     &      h_cart_norm_scale,k_cart_norm_scale)
      if (.not. active_cart_norm_scale) call errquit
     &      ('spcart_init: alloc_get failed (array scale) ',911,
     &      MEM_ERR)
      call dcopy(((2*lmaxin+1)*(lmaxin+1)),0.0d00,0,
     &    dbl_mb(k_cart_norm_scale),1)
* set normalize scale vector so cartesian vectors of transformation 
* matrix are normalized to unity
      call xlm_cart_norm_scale(lmaxin,dbl_mb(k_sp2c),
     &      dbl_mb(k_cart_norm_scale))
* set up pointers and copy recursion array to compressed 
* transformation arrays
      call xlm_ptrs(lmaxin,dbl_mb(k_sp2c),
     &      dbl_mb(k_sp2c_cmp),l_block_size,int_mb(k_sp2c_lindx))
      if (lmaxin.ge.1) then
        call xlm_ptrs_fix_p(dbl_mb(int_mb((k_sp2c_lindx+1))),3,1)
      endif
*      do lval=0,Lmaxin
*        write(6,*)' before '
*        call spcart_print_dtrans(lval)
*      enddo
*
* form inverse array 
*
      call xlm_ptrsinv(lmaxin,
     &      dbl_mb(k_invsp2c_cmp),
     &      l_block_size,
     &      int_mb(k_invsp2c_lindx))
c
* free up recursion copy of transformation matricies
      if (.not.ma_free_heap(h_sp2c)) call errquit
     &    ('spcart_init: free heap failed (array) ',911, MEM_ERR)
      active_sp2c = .false.
      k_sp2c = 0
      h_sp2c = 0
c
*      do lval=0,Lmaxin
*        write(6,*)' after '
*        call spcart_print_dtrans(lval)
*        call spcart_print_invdtrans(lval)
*      enddo
c
      sph_cart_init = Sph_Cart_Init_Value
      lmax_init = lmaxin
      spcart_init = .true.
c
Cedo#if defined(LINUX)
Cedo      trust_dgemm = .true.
Cedo#endif
c
      end
*.......................................................................
      logical function spcart_terminate()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
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
        if (active_invsp2c_cmp) then
          spcart_terminate = spcart_terminate .and.
     &        ma_free_heap(h_invsp2c_cmp)
          active_invsp2c_cmp = .false.
          k_invsp2c_cmp = 0
          h_invsp2c_cmp = 0
        endif
        if (active_invsp2c_lindx) then
          spcart_terminate = spcart_terminate .and.
     &        ma_free_heap(h_invsp2c_lindx)
          active_invsp2c_lindx = .false.
          k_invsp2c_lindx = 0
          h_invsp2c_lindx = 0
        endif
        if (active_cart_norm_scale) then
          spcart_terminate = spcart_terminate .and.
     &        ma_free_heap(h_cart_norm_scale)
          active_cart_norm_scale = .false.
          k_cart_norm_scale = 0
          h_cart_norm_scale = 0
        endif
        if (.not.spcart_terminate) call errquit
     &      (' error freeing heap in spcart_terminate',911, MEM_ERR)
        sph_cart_init = 0
        trust_dgemm = .false.
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
      subroutine xlm_ptrs(Ld,D,Dcmp,ldcmp,Dindex)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "spcartP.fh"
c
c::passed
      integer Ld        ! [input] Lmax for how spcart_* was initialized
      integer ldcmp     ! [input] length of compressed transformation 
*                                 arrays
*
*. . . . . . . . . . . . . . . . . . . . ! [input] transformation matrix
      double precision D((((Ld+1)*(Ld+2))/2),-Ld:Ld,0:Ld) 
c
      double precision Dcmp(ldcmp)  ! [output] compressed transformation
*                                              arrays
*
      integer Dindex(0:Ld)          ! [output] pointer for lth transform
*                                              array in compressed array
c
c::local      
      integer lval, isp, icart
      integer icount
c
      if (sph_cart_allsph) call errquit
     &      ('xlm_ptrs: all spherical components not working yet',911,
     &      MEM_ERR)
c
c
      icount = 0
      do 00100 lval = 0,Ld
        Dindex(lval) = k_sp2c_cmp + icount ! set pointer in index array
        do 00200 isp = -lval,lval
          do 00300 icart = 1,(((lval+1)*(lval+2))/2)
            icount = icount + 1
            Dcmp(icount) = D(icart,isp,lval) ! form separate D(xyz,sph)
*                                              arrays
00300     continue
00200   continue
00100 continue
*      call spcart_print_both(D,ld)
      end
*.......................................................................
* set up pointer information for inverse array
      subroutine xlm_ptrsinv(Ld,Dinvcmp,lDinvcmp,Dindex)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "spcartP.fh"
c
c::passed
      integer Ld       ! [input] Lmax for how spcart_* was initialized
      integer lDinvcmp ! [input] length of compressed inverse 
*                                transformation arrays
*
      double precision Dinvcmp(lDinvcmp) ! [output] compressed inverse 
*                                                   transformation 
*                                                   arrays
*
      integer Dindex(0:Ld) ! [output] pointer for lth inverse transform
*                                     array in compressed array
c
c::local      
      integer lval, ld2s, icount, l2s
      integer k_ov, h_ov, k_scr, h_scr
c
      if (sph_cart_allsph) call errquit
     &      ('xlm_ptrs: all spherical components not working yet',911,
     &      UNKNOWN_ERR)
c
      ld2s = (((ld+1)*(ld+2))/2)
      if (.not.ma_push_get(mt_dbl,ld2s*ld2s,
     &    'xlm_ptrsinv overlap',h_ov,k_ov)) call errquit
     &    ('xlm_ptrsinv: could not allocate overlap',911, MEM_ERR)
      call dcopy((ld2s*ld2s),0.0d00,0,dbl_mb(k_ov),1)
      if (.not.ma_push_get(mt_dbl,ld2s*(2*ld+1),
     &    'xlm_ptrsinv scratch',h_scr,k_scr)) call errquit
     &    ('xlm_ptrsinv: could not allocate scratch',911, MEM_ERR)
      call dcopy((ld2s*(2*ld+1)),0.0d00,0,dbl_mb(k_scr),1)
      icount = 1
      do 00100 lval = 0,Ld
        Dindex(lval) =
     &        k_invsp2c_cmp + icount - 1 ! set pointer in index array
        l2s = (((lval+1)*(lval+2))/2)
        call spcart_cart_overlap(lval,l2s,dbl_mb(k_ov))
*        write(6,*)' xlm_ptrsinv verify 1', l2s, lval
*        if (.not. ma_verify_allocator_stuff()) stop ' xlm_ptrsinv'
        call xlm_ptrsinva(l2s,lval,dbl_mb(k_ov),Dinvcmp(icount),
     &      dbl_mb(k_scr))
*        write(6,*)' xlm_ptrsinv verify 2', l2s, lval
*        if (.not. ma_verify_allocator_stuff()) stop ' xlm_ptrsinv'
        icount = icount + l2s*(2*lval+1)
00100 continue
      if (.not.ma_pop_stack(h_scr)) call errquit
     &      ('xlm_ptrsinv: could not pop_stack for overlap',911,
     &      MEM_ERR)
      if (.not.ma_pop_stack(h_ov)) call errquit
     &      ('xlm_ptrsinv: could not pop_stack for overlap',911,
     &      MEM_ERR)
      end
*.......................................................................
      subroutine xlm_ptrsinva(l2s,lval,overlap,Dinv,scr)
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"
      integer l2s,lval
      double precision Dinv(-lval:lval,l2s)
      double precision scr(-lval:lval,l2s)
      double precision overlap(l2s,l2s)
      integer s,c1, c2
      double precision val
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &      dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &      ((iis+iil)*(iil+1)*(iil+2)/2)
     &      + iic - 1)
c::statement function ----- end

* transpose Dtrans
      do s = -lval,lval
        do c1 = 1,l2s
          scr(s,c1) = Dtrans(c1,s,lval)
          Dinv(s,c1) = 0.0d00
        enddo
      enddo
*      write(6,*)' inside forming Dinv '
*      call spcart_print_dtrans(lval)
*      write(6,*)' transpose scr '
*      call output(scr,1,(2*lval+1),1,l2s,(2*lval+1),l2s,1)
* Dinv(sph,cart) = [Dtrans(cart,sph]^+ overlap(cart,cart)
* Dinv(sph,cart) = scr(sph,cart)*overlap(cart,cart)
      do s = -lval,lval
        do c1 = 1,l2s
          val = 0.0d00
          do c2 = 1,l2s
             val = val + scr(s,c2)*overlap(c2,c1)
          enddo
          Dinv(s,c1) = val
        enddo
      enddo
c
*      call dgemm('n','n',
*     &    (2*lval+1),(2*lval+1),l2s,
*     &    1.0d00,
*     &    Dinv,(2*lval+1),
*     &    dbl_mb(((int_mb(k_sp2c_lindx+lval)))),l2s,
*     &    0.0d00,scr,(2*lval+1))
*      write(6,*)' dinv*d'
*      call output(scr,1,(2*lval+1),1,(2*lval+1),(2*lval+1),(2*lval+1),1)
      end
*.......................................................................
      subroutine xlm_cart_norm_scale(Ld,D,cart_scale)
      implicit none
c::passed
      integer Ld        ! [input] Lmax for how spcart_* was initialized
      double precision
     &      D((((Ld+1)*(Ld+2))/2),-Ld:Ld,0:Ld) ! [input] transformation
*                                                        matrix
      double precision cart_scale(-Ld:Ld,0:Ld)
c::local
      integer l,c,s,l2s
      double precision norm
      call dcopy ((2*Ld+1)*(Ld+1),0.0d00,0,cart_scale,1)

      do l = 0,Ld
        l2s = (l+1)*(l+2)/2
        do s = (-l),l
          norm = 0.0d00
          do c = 1,l2s
            norm = norm + D(c,s,l)*D(c,s,l)
          enddo
*          write(6,*)'norm', norm, '  l,s,c',l,s,c
*          if (norm.gt.1.0d-10) then
          norm = sqrt(1.0d00/norm)
*          else
*            write(6,*)'scarry norm', norm, '  l,s,c',l,s,c
*            norm = 1.0d00
*          endif
          cart_scale(s,l) = norm
        enddo
      enddo
      end
*.......................................................................
*rak:      subroutine xlm_ptrs_phase(Dp,l2p,lp)
*rak:      implicit none
*rak:      integer l2p, lp
*rak:      double precision Dp(l2p,-lp:lp)
*rak:      integer lc, ls
*rak:      logical scale_it
*rak:      double precision dmaxval
*rak:      integer dmaxindx
*rak:c
*rak:      do ls = -lp,lp
*rak:        scale_it = .false.
*rak:        dmaxval  = abs(Dp(1,ls))
*rak:        dmaxindx = 1
*rak:        do lc = 2,l2p
*rak:          if (dmaxval.lt.abs(Dp(lc,ls))) then
*rak:            dmaxval = abs(Dp(lc,ls))
*rak:            dmaxindx = lc
*rak:          endif
*rak:        enddo
*rak:        if (Dp(dmaxindx,ls).lt.0.0d00) scale_it = .true.
*rak:c
*rak:        if (scale_it) then
*rak:          do lc = 1,l2p
*rak:            Dp(lc,ls) = -1.0d00*Dp(lc,ls)
*rak:          enddo
*rak:        endif
*rak:      enddo
*rak:      end
*.......................................................................
      subroutine xlm_ptrs_fix_p(Dp,l2p,lp)
      implicit none
#include "errquit.fh"
      integer l2p,lp
      integer count_val, lc, ls
      double precision Dp(l2p,-lp:lp)
c
      double precision dpdp(3)
c
      count_val = 0
      do ls = -lp,lp
        do lc = 1,l2p
          if (Dp(lc,ls).ne.0.0d00) then
            count_val = count_val + 1
            if (count_val.gt.3) then
              write(6,*)' count_val range is 1<3' 
              write(6,*)' count_val out of range ',count_val
              call errquit('fix p: error',911, MEM_ERR)
            endif
            dpdp(count_val) = Dp(lc,ls)
*            if (dpdp(count_val).lt.0.0d00)
*     &            dpdp(count_val) = -dpdp(count_val)
          endif
          Dp(lc,ls) = 0.0d00
        enddo
      enddo
      Dp(1,-1) = dpdp(1)
      Dp(2, 0) = dpdp(2)
      Dp(3, 1) = dpdp(3)
c
      end
*.......................................................................
*rak:      subroutine spcart_print_both(D,ld)
*rak:      implicit none
*rak:#include "mafdecls.fh"
*rak:#include "errquit.fh"
*rak:#include "spcartP.fh"
*rak:      integer ld
*rak:      double precision D((((Ld+1)*(Ld+2))/2),-Ld:Ld,0:Ld) ! [input] transformation matrix
*rak:c
*rak:      integer lval,l2s, ic, is
*rak:      double precision diff
*rak:c::statement function ----- start
*rak:      integer iic,iis,iil
*rak:      double precision Dtrans
*rak:      Dtrans(iic,iis,iil) =
*rak:     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
*rak:     &           ((iis+iil)*(iil+1)*(iil+2)/2)
*rak:     &           + iic - 1)
*rak:c::statement function ----- end
*rak:c
*rak:      do lval = 0,ld
*rak:        write(6,*)' d matrix '
*rak:        l2s = (lval+1)*(lval+2)/2
*rak:        do is = -lval,lval
*rak:          do ic = 1,l2s
*rak:            diff = D(ic,is,lval)-Dtrans(ic,is,lval)
*rak:            write(6,*) 
*rak:     &            lval,is,ic,D(ic,is,lval),Dtrans(ic,is,lval),diff
*rak:          enddo
*rak:        enddo
*rak:      enddo
*rak:      end
*.......................................................................
      subroutine spcart_print_dtrans(ld)
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"
      integer ld
c
      write(6,*)' spcart trans matrix used for Lval =',ld
      call output(dbl_mb((int_mb(k_sp2c_lindx+ld))),1,
     &    ((ld+1)*(ld+2)/2),1,(2*ld+1),
     &    ((ld+1)*(ld+2)/2),(2*ld+1),1)
      end
*.......................................................................
      subroutine spcart_print_invdtrans(ld)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
      integer ld
c
      write(6,*)' spcart inverse trans matrix used for Lval =',ld
      call output(dbl_mb((int_mb(k_invsp2c_lindx+ld))),1,
     &    (2*ld+1),1,((ld+1)*(ld+2)/2),
     &    (2*ld+1),((ld+1)*(ld+2)/2),1)
      end
*.......................................................................
      integer function spcart_trns_ptr(lval,lcart,lsp)
c return pointer in ma to transformation matrix for lval
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "spcartP.fh"

      integer lval  ! [input] l-value of requested matrix
      integer lcart ! [output] number of cartesian components
      integer lsp   ! [output] number of spherical components
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_trns_ptr: spcart not initialized properly',
     &    sph_cart_init, MEM_ERR)
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

      data h_sp2c            /0/
      data k_sp2c            /0/
      data h_sp2c_cmp        /0/
      data k_sp2c_cmp        /0/
      data h_invsp2c_cmp     /0/
      data k_invsp2c_cmp     /0/
      data lmax_init         /0/
      data h_cart_norm_scale /0/
      data k_cart_norm_scale /0/
      data sph_cart_init     /0/
      data h_sp2c_lindx      /0/
      data k_sp2c_lindx      /0/
      data h_invsp2c_lindx   /0/
      data k_invsp2c_lindx   /0/
      data sph_cart_allsph        /.false./
      data active_sp2c            /.false./
      data active_sp2c_cmp        /.false./
      data active_invsp2c_cmp     /.false./
      data active_sp2c_lindx      /.false./
      data active_invsp2c_lindx   /.false./
      data active_cart_norm_scale /.false./
      data trust_dgemm            /.false./
      end
*.......................................................................
      subroutine spcart_a_s(blockin, blockout, ndima, ls,
     &      ngls, in_place, print_info)
      implicit none
c
c  transforms a block of integrals with the Ls function is the slowest 
c  dimension: e.g., blockin(ndima,L2s,ngls)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ls     ! [input] angular momentum of block
      integer ngls   ! [input] general contraction length of ls block
      double precision
     &      blockin (ndima,((ls+1)*(ls+2)/2),ngls) ! [input] matrix 
      double precision
     &      blockout(ndima,-ls:ls,ngls) ! [output] matrix 
      logical
     &      in_place ! [input] true if blockin and block out are 
*                              the same pointer
      logical print_info ! [input] print info
c: local
      double precision sumg
      integer i,j,k,g
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
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (ls.lt.2) then
c...        ((ls+1)*(ls+2)/2)) = 2*ls + 1
        call dcopy((ndima*(2*ls+1)*ngls),blockin,1,blockout,1)
        return
CBERT else
CBERT   call dcopy((ndima*(2*ls+1)*ngls),0.0d00,0,blockout,1)
      endif
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
*      call spcart_print_Dtrans(ls)
c
      if (ls.eq.2) then
        call spcart_a_sd(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.3) then
        call spcart_a_sf(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.4) then
        call spcart_a_sg(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.5) then
        call spcart_a_sh(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      endif
c
c
c:old<dgemm call for blockin(ndima,l2s)*dtrans(l2s,2l+1) = 
*                                                  blockout(ndima,2l+1)>
c dgemm call for blockin(ndima,l2s,ngls)*dtrans(l2s,2l+1) = 
*                                             blockout(ndima,2l+1,ngls)
c
      L2s = ((ls+1)*(ls+2)/2)
      if (trust_dgemm) then
        do  g=1,ngls
          call dgemm('n','n',ndima,(2*ls+1),L2s,1.0d00,
     &          blockin(1,1,g),ndima,
     &          dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,
     &          0.0d00,blockout(1,1,g),ndima)
        enddo
        return
      endif
c
c
c  blockout(ndima,2l+1,ngls) = 
*                       blockin(ndima,l2s,ngls)*d_spcart(l2s,2l+1,lp=ls)
c      
      call dcopy((ndima*(2*ls+1)*ngls),0.0d00,0,blockout,1)
      L2s = ((ls+1)*(ls+2)/2)
      do g = 1,ngls
        do j=(-ls),ls
          do i=1,ndima
            sumg = 0.0d00
            do k = 1,L2s
              sumg = sumg + blockin(i,k,g)*Dtrans(k,j,Ls)
            enddo
            blockout(i,j,g) = sumg
          enddo
        enddo
      enddo
      end
*.......................................................................
C>
C> \brief Transform a block of 1-electron integrals with angular 
C> momentum <i>ls</i>
C>
C> This routine transforms integrals of a specific angular momentum 
C> <i>ls</i>, where <i>ls</i> specifies the leading dimension. As the
C> required transformation is different for different angular momenta
C> this routine essentially figures out which transformation routine to
C> call for the given angular momentum.
C>
      subroutine spcart_s_a(blockin, blockout, ndima,
     &      ls, ngls, in_place, print_info)
      implicit none
c
c  transforms a block of "integrals" with the ls function is the 
c  leading dimension; e.g., blockin(L2s,ndima)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  !< [Input] leading dimension of block
      integer ls     !< [Input] angular momentum of block
      integer ngls   !< [Input] general contraction length of ls block
      double precision
     &      blockin (((ls+1)*(ls+2)/2),ngls,ndima) !< [Input] Cartesian
                                                   !< integrals
      double precision
     &      blockout(-ls:ls,ngls,ndima) !< [Output] spherical harmonic
                                        !< integrals 
      logical in_place   !< [Input] true if blockin and block out are 
                         !< the same pointer
      logical print_info !< [Input] print info
c: local
      integer i,j,k,g
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
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (ls.lt.2) then
c...        ((ls+1)*(ls+2)/2)) = 2*ls + 1
        call dcopy((ndima*(2*ls+1)*ngls),blockin,1,blockout,1)
        return
CBERT else
CBERT   call dcopy((ndima*(2*ls+1)*ngls),0.0d00,0,blockout,1)
      endif
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
*      call spcart_print_Dtrans(ls)
c
      if (ls.eq.2) then
        call spcart_sd_a(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.3) then
        call spcart_sf_a(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.4) then
        call spcart_sg_a(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.5) then
        call spcart_sh_a(blockin,blockout,ndima,ngls,
     &      in_place,print_info)
        return
      endif
      if (ngls.eq.1 .and. trust_dgemm) then
c
c only works for ngls = 1 right now
c dgemm call for:
c  Transpose(d_spcart(l2s,2l+1))*blockin(l2s,ndima) = 
*                                                   blockout(2l+1,ndima)
c
        L2s = ((ls+1)*(ls+2)/2)
        call dgemm('t','n',(2*ls+1),ndima,L2s,1.0d00,
     &        dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,
     &        blockin,L2s,0.0d00,blockout,(2*ls+1))
        return
      endif
c
c blockout(2l+1,ngls,ndima) = 
*                       blockin(l2s,ngls,ndima)*d_spcart(l2s,2l+1,lp=ls)
c      
      L2s = ((ls+1)*(ls+2)/2)
      call dcopy((ndima*(2*ls+1)*ngls),0.0d00,0,blockout,1)
      do j=1,ndima
        do g=1,ngls
          do i=(-ls),ls
            do k = 1,L2s
              blockout(i,g,j) = blockout(i,g,j) +
     &              blockin(k,g,j)*Dtrans(k,i,Ls)
            enddo
          enddo
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_a_s_b(blockin, blockout, ndima, ndimb,
     &      ls, ngls,
     &      in_place,print_info)
      implicit none
c
c  transforms a block of "integrals" with the ls function is ordered 
c  between a leading dimension and trailing dimension.
c  e.g., blockin(nidima,L2s,ndimb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] tailing dimension of block
      integer ls     ! [input] angular momentum of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,((ls+1)*(ls+2)/2),ngls,ndimb) 
*                                                      ! [input] matrix
*
      double precision
     &      blockout(ndima,-ls:ls,ngls,ndimb) ! [output] matrix 
      logical
     &      in_place     ! [input] true if blockin and blockout are 
*                                  the same pointer
*
      logical print_info ! [input] print info
c: local
      integer i,j,k,g
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
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (ls.lt.2) then
c...        ((ls+1)*(ls+2)/2)) = 2*ls + 1
        call dcopy((ndima*ndimb*(2*ls+1)*ngls),blockin,1,blockout,1)
        return
CBERT else
CBERT   call dcopy((ndima*ndimb*(2*ls+1)*ngls),0.0d00,0,blockout,1)
      endif
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
*      call spcart_print_Dtrans(ls)
c
      if (ls.eq.2) then
        call spcart_a_sd_b(blockin,blockout,ndima,ndimb,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.3) then
        call spcart_a_sf_b(blockin,blockout,ndima,ndimb,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.4) then
        call spcart_a_sg_b(blockin,blockout,ndima,ndimb,ngls,
     &      in_place,print_info)
        return
      else if (ls.eq.5) then
        call spcart_a_sh_b(blockin,blockout,ndima,ndimb,ngls,
     &      in_place,print_info)
        return
      endif
c
      if (ngls.eq.1 .and. trust_dgemm) then
c
c dgemm for all ndimb only if ngls == 1
c
        L2s = ((ls+1)*(ls+2)/2)
        do m = 1,ndimb
          call dgemm('n','n',ndima,(2*ls+1),L2s,
     &          1.0d00,blockin(1,1,1,m),ndima,
     &          dbl_mb(int_mb(k_sp2c_lindx+Ls)),L2s,
     &          0.0d00,blockout(1,1,1,m),ndima)
        enddo
        return
      endif
c
      call dcopy((ndima*ndimb*(2*ls+1)*ngls),0.0d00,0,blockout,1)
      L2s = ((ls+1)*(ls+2)/2)
      do  m = 1,ndimb
        do  g=1,ngls
          do  j=(-ls),ls
            do  k = 1,L2s
              do  i=1,ndima
                blockout(i,j,g,m) = blockout(i,j,g,m) +
     &                blockin(i,k,g,m)*Dtrans(k,j,Ls)
              enddo
            enddo
          enddo
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_a_sd(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the d function is the slowest 
c  dimension; e.g., blockin(ndima,6) and blockout(ndima,5)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,6,ngls)    ! [input] matrix 
      double precision blockout(ndima,-2:2,ngls) ! [output] matrix 
      logical
     &      in_place     ! [input] true if blockin and block out are 
*                                  the same pointer
*
      logical print_info ! [input] print info
c: local
      integer i,g
*rak:      integer ls     ! [fixed at 2] angular momentum of block
      double precision dt2m2, dt5m1, dt10, dt40, dt60, dt31, dt12, dt42
      double precision bin1, bin2, bin3, bin4, bin5, bin6
      double precision sinm2, sinm1, sin0, sinp1, sinp2
      logical print_debug
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
      print_debug = .true.
c
      if (print_info.and.print_debug) then
        write(6,*)' insize spcart_a_sd'
        write(6,*)' spcart_a_sd:ndima = ',ndima
        write(6,*)' spcart_a_sd:ngls  = ',ngls
      endif
      if (print_info) call spcart_print_dtrans(2)
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
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c  spint(-2) = cint(2)*dtrans(2,-2,2)
c  spint(-1) = cint(5)*dtrans(5,-1,2)
c  spint( 0) = cint(1)*dtrans(1,0,2) + cint(4)*dtrans(4,0,2) 
*              + cint(6)*dtrans(6,0,2)
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
      do g=1,ngls
        do i=1,ndima
#ifdef DEBUG
          if (print_info.and.print_debug)
     &        write(6,*)' spcart_a_sd:g,i',g,i
#endif
          bin1 = blockin(i,1,g)
          bin4 = blockin(i,4,g)
#ifdef DEBUG
          if (print_info.and.print_debug) then   
            write(6,11111)' spcart_a_sd:cart ints 1 xx',bin1,g,i
            write(6,11111)' spcart_a_sd:cart ints 2 xy',bin2,g,i
            write(6,11111)' spcart_a_sd:cart ints 3 xz',bin3,g,i
            write(6,11111)' spcart_a_sd:cart ints 4 yy',bin4,g,i
            write(6,11111)' spcart_a_sd:cart ints 6 zz',bin6,g,i
          endif
#endif
          blockout(i,-2,g) = blockin(i,2,g)*dt2m2
          blockout(i,-1,g) = blockin(i,5,g)*dt5m1
          blockout(i, 0,g) = bin1*dt10+bin4*dt40+blockin(i,6,g)*dt60
          blockout(i, 1,g) = blockin(i,3,g)*dt31
          blockout(i, 2,g) = bin1*dt12+bin4*dt42
#ifdef DEBUG
          if (print_info.and.print_debug) then   
            write(6,11111)' spcart_a_sd:sph ints -2',sinm2,g,i
            write(6,11111)' spcart_a_sd:sph ints -1',sinm1,g,i
            write(6,11111)' spcart_a_sd:sph ints  0',sin0,g,i
            write(6,11111)' spcart_a_sd:sph ints +1',sinp1,g,i
            write(6,11111)' spcart_a_sd:sph ints +2',sinp2,g,i
          endif
#endif
        enddo
      enddo
      if (print_info.and.print_debug)
     &    write(6,*)' exiting spcart_a_sd'
11111 format(1x,a,1pd20.10,i5,i5)
      end
*.......................................................................
C>
C> \brief Transforms a block of 1-electron integrals of D-functions
C>
C> This routine transforms a block of integrals where the functions
C> corresponding to the leading dimension are Cartesian D-functions
C> on input and spherical harmonic D-functions on output.
C>
      subroutine spcart_sd_a(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the d function is the fastest 
c  dimension; e.g., blockin(6,ndima) and blockout(5,ndima)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  !< [Input] leading dimension of block
      integer ngls   !< [Input] general contraction length of ls block
      double precision blockin (6,ngls,ndima)    !< [Input] matrix 
      double precision blockout(-2:2,ngls,ndima) !< [Output] matrix 
      logical in_place   !< [Input] true if blockin and block out are 
                         !< the same pointer
      logical print_info !< [Input] print info
c: local
      integer i,g
*rak:      integer ls     ! [fixed at 2] angular momentum of block
      double precision dt2m2, dt5m1, dt10, dt40, dt60, dt31, dt12, dt42
      double precision bin1, bin2, bin3, bin4, bin5, bin6
      double precision sinm2, sinm1, sin0, sinp1, sinp2
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
      if (print_info) then
        write(6,*)' insize spcart_sd_a'
        write(6,*)' spcart_sd_a:ndima = ',ndima
        write(6,*)' spcart_sd_a:ngls  = ',ngls
        call spcart_print_dtrans(2)
      endif
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
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c  spint(-2) = cint(2)*dtrans(2,-2,2)
c  spint(-1) = cint(5)*dtrans(5,-1,2)
c  spint( 0) = cint(1)*dtrans(1,0,2) + cint(4)*dtrans(4,0,2) 
*              + cint(6)*dtrans(6,0,2)
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
      do i=1,ndima
        do g=1,ngls
          bin1 = blockin(1,g,i)
          bin2 = blockin(2,g,i)
          bin3 = blockin(3,g,i)
          bin4 = blockin(4,g,i)
          bin5 = blockin(5,g,i)
          bin6 = blockin(6,g,i)
          if (print_info) then
            write(6,11111)' spcart_sd_a:cart ints 1 xx',bin1,g,i
            write(6,11111)' spcart_sd_a:cart ints 2 xy',bin2,g,i
            write(6,11111)' spcart_sd_a:cart ints 3 xz',bin3,g,i
            write(6,11111)' spcart_sd_a:cart ints 4 yy',bin4,g,i
            write(6,11111)' spcart_sd_a:cart ints 5 yz',bin5,g,i
            write(6,11111)' spcart_sd_a:cart ints 6 zz',bin6,g,i
          endif
          sinm2 = bin2*dt2m2
          sinm1 = bin5*dt5m1
          sin0  = bin1*dt10 +
     &            bin4*dt40 +
     &            bin6*dt60
          sinp1 = bin3*dt31
          sinp2 = bin1*dt12 +
     &            bin4*dt42
          blockout(-2,g,i) = sinm2
          blockout(-1,g,i) = sinm1
          blockout( 0,g,i) = sin0
          blockout( 1,g,i) = sinp1
          blockout( 2,g,i) = sinp2
          if (print_info) then
            write(6,11111)' spcart_sd_a:sph ints -2',sinm2,g,i
            write(6,11111)' spcart_sd_a:sph ints -1',sinm1,g,i
            write(6,11111)' spcart_sd_a:sph ints  0',sin0,g,i
            write(6,11111)' spcart_sd_a:sph ints +1',sinp1,g,i
            write(6,11111)' spcart_sd_a:sph ints +2',sinp2,g,i
          endif
        enddo
      enddo
      if (print_info)
     &    write(6,*)' exiting spcart_sd_a'
11111 format(1x,a,1pd20.10,i5,i5)
      end
*.......................................................................
      subroutine spcart_a_sd_b(blockin, blockout,
     &    ndima, ndimb, ngls, in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the d function is the middle 
c  dimension; e.g., blockin(ndima,6,ndimb) and blockout(ndima,5,ndimb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] trailing dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision
     &      blockin (ndima, 6   ,ngls,ndimb) ! [input] matrix 
      double precision
     &      blockout(ndima, -2:2,ngls,ndimb) ! [output] matrix 
      logical
     &      in_place     ! [input] true if blockin and block out are 
*                                  the same pointer
      logical print_info ! [input] print info
c: local
      integer i,j,g
*rak:      integer ls     ! [fixed at 2] angular momentum of block
      double precision dt2m2, dt5m1, dt10, dt40, dt60, dt31, dt12, dt42
      double precision bin1, bin2, bin3, bin4, bin5, bin6
      double precision sinm2, sinm1, sin0, sinp1, sinp2
c::statement function ----- start
      integer iic,iis,iil,i1
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
      if (print_info) then
        write(6,*)' insize spcart_a_sd_b'
        write(6,*)' spcart_a_sd_b:ndima = ',ndima
        write(6,*)' spcart_a_sd_b:ngls  = ',ngls
        write(6,*)' spcart_a_sd_b:ndimb = ',ndimb
        call spcart_print_dtrans(2)
      endif
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
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c  spint(-2) = cint(2)*dtrans(2,-2,2)
c  spint(-1) = cint(5)*dtrans(5,-1,2)
c  spint( 0) = cint(1)*dtrans(1,0,2) + cint(4)*dtrans(4,0,2) 
*              + cint(6)*dtrans(6,0,2)
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
      do j=1,ndimb
        do g=1,ngls
#define CHUNK 4
!DEC$ LOOP COUNT AVG=CHUNK
          do i1=1,ndima,CHUNK
          do i=i1,min(i1+CHUNK-1,ndima)
#ifdef DEBUG
            if (print_info)
     &          write(6,*)' spcart_a_sd_b:j,g,i',j,g,i
#endif
            bin1 = blockin(i,1,g,j)
            bin2 = blockin(i,2,g,j)
            bin3 = blockin(i,3,g,j)
            bin4 = blockin(i,4,g,j)
            bin5 = blockin(i,5,g,j)
            bin6 = blockin(i,6,g,j)
#ifdef DEBUG
            if (print_info) then
              write(6,11111)' spcart_a_sd_b:cart ints 1 xx',bin1,j,g,i
              write(6,11111)' spcart_a_sd_b:cart ints 2 xy',bin2,j,g,i
              write(6,11111)' spcart_a_sd_b:cart ints 3 xz',bin3,j,g,i
              write(6,11111)' spcart_a_sd_b:cart ints 4 yy',bin4,j,g,i
              write(6,11111)' spcart_a_sd_b:cart ints 5 yz',bin5,j,g,i
              write(6,11111)' spcart_a_sd_b:cart ints 6 zz',bin6,j,g,i
            endif
#endif
            sinm2 = bin2*dt2m2
            sinm1 = bin5*dt5m1
            sin0  = bin1*dt10 +
     &              bin4*dt40 +
     &              bin6*dt60
            sinp1 = bin3*dt31
            sinp2 = bin1*dt12 +
     &              bin4*dt42
            blockout(i,-2,g,j) = sinm2
            blockout(i,-1,g,j) = sinm1
            blockout(i, 0,g,j) = sin0
            blockout(i, 1,g,j) = sinp1
            blockout(i, 2,g,j) = sinp2
#ifdef DEBUG
            if (print_info) then
              write(6,11111)' spcart_a_sd_b:sph ints -2',sinm2,j,g,i
              write(6,11111)' spcart_a_sd_b:sph ints -1',sinm1,j,g,i
              write(6,11111)' spcart_a_sd_b:sph ints  0',sin0,j,g,i
              write(6,11111)' spcart_a_sd_b:sph ints +1',sinp1,j,g,i
              write(6,11111)' spcart_a_sd_b:sph ints +2',sinp2,j,g,i
            endif
#endif
          enddo
          enddo
        enddo
      enddo
      if (print_info)
     &    write(6,*)' exiting spcart_a_sd_b'
11111 format(1x,a,1pd20.10,i5,i5,i5)
      end
*.......................................................................
      subroutine spcart_a_sf(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the f function is the slowest 
c  dimension; e.g., blockin(ndima,10) and blockout(ndima,7)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,10,ngls)    ! [input] matrix 
      double precision blockout(ndima,-3:3,ngls)  ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are 
*                                 the same pointer
      logical print_info ! [input] print info
c: local
      integer i,g
*rak:      integer ls     ! [fixed at 3] angular momentum of block
      double precision dt2_m3, dt7_m3, dt5_m2, dt2_m1, dt7_m1, dt9_m1
      double precision dt3_0, dt8_0, dt10_0, dt1_1, dt4_1, dt6_1
      double precision dt3_2, dt8_2, dt1_3, dt4_3
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
c::statement function ----- start
      integer iic,iis,iil,i1
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
     &    sph_cart_init, UNKNOWN_ERR)
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
      do g=1,ngls
#define CHUNK 4
!DEC$ LOOP COUNT AVG=CHUNK
          do i1=1,ndima,CHUNK
          do i=i1,min(i1+CHUNK-1,ndima)
          bin1 = blockin(i,1,g)
          bin2 = blockin(i,2,g)
          bin3 = blockin(i,3,g)
          bin4 = blockin(i,4,g)
          bin5 = blockin(i,5,g)
          bin6 = blockin(i,6,g)
          bin7 = blockin(i,7,g)
          bin8 = blockin(i,8,g)
          bin9 = blockin(i,9,g)
          bin10 = blockin(i,10,g)
          blockout(i,-3,g) = bin2*dt2_m3 + bin7*dt7_m3
          blockout(i,-2,g) = bin5*dt5_m2
          blockout(i,-1,g) = bin2*dt2_m1 + bin7*dt7_m1 + bin9*dt9_m1
          blockout(i, 0,g) = bin3*dt3_0 + bin8*dt8_0 + bin10*dt10_0
          blockout(i, 1,g) = bin1*dt1_1 + bin4*dt4_1 + bin6*dt6_1
          blockout(i, 2,g) = bin3*dt3_2 + bin8*dt8_2
          blockout(i, 3,g) = bin1*dt1_3 + bin4*dt4_3
        enddo
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_sf_a(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the f function is the fastest dimension
c  e.g., blockin(10,ndima) and blockout(7,ndima)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (10,ngls,ndima)    ! [input] matrix 
      double precision blockout(-3:3,ngls,ndima)  ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,g
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
     &    sph_cart_init, UNKNOWN_ERR)
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
      do i=1,ndima
        do g=1,ngls
          bin1  = blockin( 1,g,i)
          bin2  = blockin( 2,g,i)
          bin3  = blockin( 3,g,i)
          bin4  = blockin( 4,g,i)
          bin5  = blockin( 5,g,i)
          bin6  = blockin( 6,g,i)
          bin7  = blockin( 7,g,i)
          bin8  = blockin( 8,g,i)
          bin9  = blockin( 9,g,i)
          bin10 = blockin(10,g,i)
          blockout(-3,g,i) = bin2*dt2_m3 + bin7*dt7_m3
          blockout(-2,g,i) = bin5*dt5_m2
          blockout(-1,g,i) = bin2*dt2_m1 + bin7*dt7_m1 + bin9*dt9_m1
          blockout( 0,g,i) = bin3*dt3_0 + bin8*dt8_0 + bin10*dt10_0
          blockout( 1,g,i) = bin1*dt1_1 + bin4*dt4_1 + bin6*dt6_1
          blockout( 2,g,i) = bin3*dt3_2 + bin8*dt8_2
          blockout( 3,g,i) = bin1*dt1_3 + bin4*dt4_3
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_a_sf_b(blockin, blockout,
     &    ndima, ndimb, ngls, in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the f function is the middle dimension
c  e.g., blockin(ndima,10,ndimb) and blockout(ndima,7,ndimb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] trailing dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,10,ngls,ndimb)    ! [input] matrix 
      double precision blockout(ndima,-3:3,ngls,ndimb)  ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,j,g
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
     &    sph_cart_init, UNKNOWN_ERR)
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
      do j=1,ndimb
        do g=1,ngls
          do i=1,ndima
            bin1  = blockin(i, 1,g,j)
            bin2  = blockin(i, 2,g,j)
            bin3  = blockin(i, 3,g,j)
            bin4  = blockin(i, 4,g,j)
            bin5  = blockin(i, 5,g,j)
            bin6  = blockin(i, 6,g,j)
            bin7  = blockin(i, 7,g,j)
            bin8  = blockin(i, 8,g,j)
            bin9  = blockin(i, 9,g,j)
            bin10 = blockin(i,10,g,j)
            blockout(i,-3,g,j) = bin2*dt2_m3 + bin7*dt7_m3
            blockout(i,-2,g,j) = bin5*dt5_m2
            blockout(i,-1,g,j) =
     &                         bin2*dt2_m1 + bin7*dt7_m1 + bin9*dt9_m1
            blockout(i, 0,g,j) =
     &                         bin3*dt3_0 + bin8*dt8_0 + bin10*dt10_0
            blockout(i, 1,g,j) =
     &                         bin1*dt1_1 + bin4*dt4_1 + bin6*dt6_1
            blockout(i, 2,g,j) = bin3*dt3_2 + bin8*dt8_2
            blockout(i, 3,g,j) = bin1*dt1_3 + bin4*dt4_3
          enddo
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_a_sg(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the g function is the slowest dimension
c  e.g., blockin(ndima,15) and blockout(ndima,9)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,15,ngls)   ! [input] matrix 
      double precision blockout(ndima,-4:4,ngls) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,g
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
     &    sph_cart_init, UNKNOWN_ERR)
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
      do g=1,ngls
        do i=1,ndima
          bin1  = blockin(i, 1,g)
          bin2  = blockin(i, 2,g)
          bin3  = blockin(i, 3,g)
          bin4  = blockin(i, 4,g)
          bin5  = blockin(i, 5,g)
          bin6  = blockin(i, 6,g)
          bin7  = blockin(i, 7,g)
          bin8  = blockin(i, 8,g)
          bin9  = blockin(i, 9,g)
          bin10 = blockin(i,10,g)
          bin11 = blockin(i,11,g)
          bin12 = blockin(i,12,g)
          bin13 = blockin(i,13,g)
          bin14 = blockin(i,14,g)
          bin15 = blockin(i,15,g)
c
          blockout(i,-4,g) = bin2*dt2_m4  + bin7*dt7_m4  
          blockout(i,-3,g) = bin5*dt5_m3  + bin12*dt12_m3 
          blockout(i,-2,g) = bin2*dt2_m2  + bin7*dt7_m2   +
     &                       bin9*dt9_m2
          blockout(i,-1,g) = bin5*dt5_m1  + bin12*dt12_m1 +
     &                       bin14*dt14_m1
          blockout(i, 0,g) = bin1*dt1_0   + bin4*dt4_0    +
     &                       bin6*dt6_0   + bin11*dt11_0  +
     &                       bin13*dt13_0 + bin15*dt15_0
          blockout(i, 1,g) = bin3*dt3_1   + bin8*dt8_1    +
     &                       bin10*dt10_1
          blockout(i, 2,g) = bin1*dt1_2   + bin6*dt6_2    +
     &                       bin11*dt11_2 + bin13*dt13_2
          blockout(i, 3,g) = bin3*dt3_3   + bin8*dt8_3
          blockout(i, 4,g) = bin1*dt1_4   + bin4*dt4_4    +
     &                       bin11*dt11_4
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_sg_a(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the g function is the fastest dimension
c  e.g., blockin(15,ndima) and blockout(9,ndima)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (15,ngls,ndima)   ! [input] matrix 
      double precision blockout(-4:4,ngls,ndima) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,g
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
     &    sph_cart_init, UNKNOWN_ERR)
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
      do i=1,ndima
        do g=1,ngls
          bin1  = blockin( 1,g,i)
          bin2  = blockin( 2,g,i)
          bin3  = blockin( 3,g,i)
          bin4  = blockin( 4,g,i)
          bin5  = blockin( 5,g,i)
          bin6  = blockin( 6,g,i)
          bin7  = blockin( 7,g,i)
          bin8  = blockin( 8,g,i)
          bin9  = blockin( 9,g,i)
          bin10 = blockin(10,g,i)
          bin11 = blockin(11,g,i)
          bin12 = blockin(12,g,i)
          bin13 = blockin(13,g,i)
          bin14 = blockin(14,g,i)
          bin15 = blockin(15,g,i)
c
          blockout(-4,g,i) = bin2*dt2_m4  + bin7*dt7_m4  
          blockout(-3,g,i) = bin5*dt5_m3  + bin12*dt12_m3 
          blockout(-2,g,i) = bin2*dt2_m2  + bin7*dt7_m2   +
     &                       bin9*dt9_m2
          blockout(-1,g,i) = bin5*dt5_m1  + bin12*dt12_m1 +
     &                       bin14*dt14_m1
          blockout( 0,g,i) = bin1*dt1_0   + bin4*dt4_0    +
     &                       bin6*dt6_0   + bin11*dt11_0  +
     &                       bin13*dt13_0 + bin15*dt15_0
          blockout( 1,g,i) = bin3*dt3_1   + bin8*dt8_1    +
     &                       bin10*dt10_1
          blockout( 2,g,i) = bin1*dt1_2   + bin6*dt6_2    +
     &                       bin11*dt11_2 + bin13*dt13_2
          blockout( 3,g,i) = bin3*dt3_3   + bin8*dt8_3
          blockout( 4,g,i) = bin1*dt1_4   + bin4*dt4_4    +
     &                       bin11*dt11_4
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_a_sg_b(blockin, blockout,
     &    ndima, ndimb, ngls, in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the g function is the middle dimension
c  e.g., blockin(ndima,15,ndimb) and blockout(ndima,9,ndimb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] trailing dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,15,ngls,ndimb)   ! [input] matrix 
      double precision blockout(ndima,-4:4,ngls,ndimb) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,j,g
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
     &    sph_cart_init, UNKNOWN_ERR)
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
      do j=1,ndimb
        do g=1,ngls
          do i=1,ndima
            bin1  = blockin(i, 1,g,j)
            bin2  = blockin(i, 2,g,j)
            bin3  = blockin(i, 3,g,j)
            bin4  = blockin(i, 4,g,j)
            bin5  = blockin(i, 5,g,j)
            bin6  = blockin(i, 6,g,j)
            bin7  = blockin(i, 7,g,j)
            bin8  = blockin(i, 8,g,j)
            bin9  = blockin(i, 9,g,j)
            bin10 = blockin(i,10,g,j)
            bin11 = blockin(i,11,g,j)
            bin12 = blockin(i,12,g,j)
            bin13 = blockin(i,13,g,j)
            bin14 = blockin(i,14,g,j)
            bin15 = blockin(i,15,g,j)
c
            blockout(i,-4,g,j) = bin2*dt2_m4  + bin7*dt7_m4  
            blockout(i,-3,g,j) = bin5*dt5_m3  + bin12*dt12_m3 
            blockout(i,-2,g,j) = bin2*dt2_m2  + bin7*dt7_m2   +
     &                           bin9*dt9_m2
            blockout(i,-1,g,j) = bin5*dt5_m1  + bin12*dt12_m1 +
     &                           bin14*dt14_m1
            blockout(i, 0,g,j) = bin1*dt1_0   + bin4*dt4_0    +
     &                           bin6*dt6_0   + bin11*dt11_0  +
     &                           bin13*dt13_0 + bin15*dt15_0
            blockout(i, 1,g,j) = bin3*dt3_1   + bin8*dt8_1    +
     &                           bin10*dt10_1
            blockout(i, 2,g,j) = bin1*dt1_2   + bin6*dt6_2    +
     &                           bin11*dt11_2 + bin13*dt13_2
            blockout(i, 3,g,j) = bin3*dt3_3   + bin8*dt8_3
            blockout(i, 4,g,j) = bin1*dt1_4   + bin4*dt4_4    +
     &                           bin11*dt11_4
          enddo
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_a_sh(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the h function is the fastest dimension
c  e.g., blockin(ndima,21) and blockout(ndima,11)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,21,ngls)   ! [input] matrix 
      double precision blockout(ndima,-5:5,ngls) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,g
      double precision dt2_m5, dt7_m5, dt16_m5, dt5_m4, dt12_m4, dt2_m3
      double precision dt7_m3, dt9_m3, dt16_m3, dt18_m3, dt5_m2, dt12_m2
      double precision dt14_m2, dt2_m1, dt7_m1, dt9_m1, dt16_m1, dt18_m1
      double precision dt20_m1, dt3_0, dt8_0, dt10_0, dt17_0, dt19_0
      double precision dt21_0, dt1_1, dt4_1, dt6_1, dt11_1, dt13_1 
      double precision dt15_1, dt3_2, dt10_2, dt17_2, dt19_2, dt1_3
      double precision dt4_3, dt6_3, dt11_3, dt13_3, dt3_4, dt8_4
      double precision dt17_4, dt1_5, dt4_5, dt11_5
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
      double precision bin11, bin12, bin13, bin14, bin15
      double precision bin16, bin17, bin18, bin19, bin20, bin21
*rak:      integer ls     ! [fixed at 5] angular momentum of block
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=5
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_sh_a: spcart not initialized properly',
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
c  spint(-5) = cint(2)*dtrans(2,-5,5) + cint(7)*dtrans(7,-5,5) + cint(16)*dtrans(16,-5,5)
c  spint(-4) = cint(5)*dtrans(5,-4,5) + cint(12)*dtrans(12,-4,5)
c  spint(-3) = cint(2)*dtrans(2,-3,5) + cint(7)*dtrans(7,-3,5) + cint(9)*dtrans(9,-3,5)
c            + cint(16)*dtrans(16,-3,5) + cint(18)*dtrans(18,-3,5)
c  spint(-2) = cint(5)*dtrans(5,-2,5) + cint(12)*dtrans(12,-2,5) + cint(14)*dtrans(14,-2,5)
c  spint(-1) = cint(2)*dtrans(2,-1,5) + cint(7)*dtrans(7,-1,5) + cint(9)*dtrans(9,-1,5)
c            + cint(16)*dtrans(16,-1,5) + cint(18)*dtrans(18,-1,5) + cint(20)*dtrans(20,-1,5) 
c  spint( 0) = cint(3)*dtrans(3,0,5) + cint(8)*dtrans(8,0,5) + cint(10)*dtrans(10,0,5)
c            + cint(17)*dtrans(17,0,5) + cint(19)*dtrans(19,0,5) + cint(21)*dtrans(21,0,5) 
c  spint( 1) = cint(1)*dtrans(1,1,5)   + cint(4)*dtrans(4,1,5)   + cint(6)*dtrans(6,1,5)
c            + cint(11)*dtrans(11,1,5) + cint(13)*dtrans(13,1,5) + cint(15)*dtrans(15,1,5)
c  spint( 2) = cint(3)*dtrans(3,2,5) + cint(10)*dtrans(10,2,5) + cint(17)*dtrans(17,2,5)
c            + cint(19)*dtrans(19,2,5)
c  spint( 3) = cint(1)*dtrans(1,3,5) + cint(4)*dtrans(4,3,5) + cint(6)*dtrans(6,3,5)
c            + cint(11)*dtrans(11,3,5) + cint(13)*dtrans(13,3,5)
c  spint( 4) = cint(3)*dtrans(3,4,5) + cint(8)*dtrans(8,4,5) + cint(17)*dtrans(17,4,5)
c  spint( 5) = cint(1)*dtrans(1,4,5) + cint(4)*dtrans(4,4,5) + cint(11)*dtrans(11,4,5)
c
      dt2_m5  = dtrans(2,-5,5)
      dt7_m5  = dtrans(7,-5,5)
      dt16_m5 = dtrans(16,-5,5)
      dt5_m4  = dtrans(5,-4,5)     
      dt12_m4 = dtrans(12,-4,5)    
      dt2_m3  = dtrans(2,-3,5)     
      dt7_m3  = dtrans(7,-3,5)     
      dt9_m3  = dtrans(9,-3,5)     
      dt16_m3 = dtrans(16,-3,5)
      dt18_m3 = dtrans(18,-3,5)
      dt5_m2  = dtrans(5,-2,5)     
      dt12_m2 = dtrans(12,-2,5)    
      dt14_m2 = dtrans(14,-2,5)    
      dt2_m1  = dtrans(2,-1,5)
      dt7_m1  = dtrans(7,-1,5)
      dt9_m1  = dtrans(9,-1,5)
      dt16_m1 = dtrans(16,-1,5)
      dt18_m1 = dtrans(18,-1,5)
      dt20_m1 = dtrans(20,-1,5)
      dt3_0   = dtrans(3,0,5)      
      dt8_0   = dtrans(8,0,5)      
      dt10_0  = dtrans(10,0,5)     
      dt17_0  = dtrans(17,0,5)     
      dt19_0  = dtrans(19,0,5)     
      dt21_0  = dtrans(21,0,5)     
      dt1_1   = dtrans(1,1,5)      
      dt4_1   = dtrans(4,1,5)      
      dt6_1   = dtrans(6,1,5)      
      dt11_1  = dtrans(11,1,5)     
      dt13_1  = dtrans(13,1,5)     
      dt15_1  = dtrans(15,1,5)     
      dt3_2   = dtrans(3,2,5)      
      dt10_2  = dtrans(10,2,5)     
      dt17_2  = dtrans(17,2,5)     
      dt19_2  = dtrans(19,2,5)     
      dt1_3   = dtrans(1,3,5)      
      dt4_3   = dtrans(4,3,5)      
      dt6_3   = dtrans(6,3,5)      
      dt11_3  = dtrans(11,3,5)     
      dt13_3  = dtrans(13,3,5)                 
      dt3_4   = dtrans(3,4,5)      
      dt8_4   = dtrans(8,4,5)      
      dt17_4  = dtrans(17,4,5)      
      dt1_5   = dtrans(1,5,5)      
      dt4_5   = dtrans(4,5,5)      
      dt11_5  = dtrans(11,5,5)     
      do g=1,ngls
        do i=1,ndima
          bin1  = blockin(i, 1,g)
          bin2  = blockin(i, 2,g)
          bin3  = blockin(i, 3,g)
          bin4  = blockin(i, 4,g)
          bin5  = blockin(i, 5,g)
          bin6  = blockin(i, 6,g)
          bin7  = blockin(i, 7,g)
          bin8  = blockin(i, 8,g)
          bin9  = blockin(i, 9,g)
          bin10 = blockin(i,10,g)
          bin11 = blockin(i,11,g)
          bin12 = blockin(i,12,g)
          bin13 = blockin(i,13,g)
          bin14 = blockin(i,14,g)
          bin15 = blockin(i,15,g)
          bin16 = blockin(i,16,g)
          bin17 = blockin(i,17,g)
          bin18 = blockin(i,18,g)
          bin19 = blockin(i,19,g)
          bin20 = blockin(i,20,g)
          bin21 = blockin(i,21,g)
c
          blockout(i,-5,g) = bin2*dt2_m5  + bin7*dt7_m5   +
     &                       bin16*dt16_m5
          blockout(i,-4,g) = bin5*dt5_m4  + bin12*dt12_m4 
          blockout(i,-3,g) = bin2*dt2_m3  + bin7*dt7_m3   +
     &                       bin9*dt9_m3  + bin16*dt16_m3 +
     &                       bin18*dt18_m3
          blockout(i,-2,g) = bin5*dt5_m2  + bin12*dt12_m2 +
     &                       bin14*dt14_m2
          blockout(i,-1,g) = bin2*dt2_m1  + bin7*dt7_m1   +
     &                       bin9*dt9_m1  + bin16*dt16_m1 +
     &                       bin18*dt18_m1+ bin20*dt20_m1
          blockout(i, 0,g) = bin3*dt3_0   + bin8*dt8_0    +
     &                       bin10*dt10_0 + bin17*dt17_0  +
     &                       bin19*dt19_0 + bin21*dt21_0
          blockout(i, 1,g) = bin1*dt1_1   + bin4*dt4_1    +
     &                       bin6*dt6_1   + bin11*dt11_1  +
     &                       bin13*dt13_1 + bin15*dt15_1
          blockout(i, 2,g) = bin3*dt3_2   + bin10*dt10_2  +
     &                       bin17*dt17_2 + bin19*dt19_2
          blockout(i, 3,g) = bin1*dt1_3   + bin4*dt4_3    +
     &                       bin6*dt6_3   + bin11*dt11_3  +
     &                       bin13*dt13_3
          blockout(i, 4,g) = bin3*dt3_4   + bin8*dt8_4    +
     &                       bin17*dt17_4
          blockout(i, 5,g) = bin1*dt1_5   + bin4*dt4_5    +
     &                       bin11*dt11_5
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_sh_a(blockin, blockout, ndima, ngls,
     &    in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the h function is the fastest dimension
c  e.g., blockin(21,ndima) and blockout(11,ndima)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (21,ngls,ndima)   ! [input] matrix 
      double precision blockout(-5:5,ngls,ndima) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,g
      double precision dt2_m5, dt7_m5, dt16_m5, dt5_m4, dt12_m4, dt2_m3
      double precision dt7_m3, dt9_m3, dt16_m3, dt18_m3, dt5_m2, dt12_m2
      double precision dt14_m2, dt2_m1, dt7_m1, dt9_m1, dt16_m1, dt18_m1
      double precision dt20_m1, dt3_0, dt8_0, dt10_0, dt17_0, dt19_0
      double precision dt21_0, dt1_1, dt4_1, dt6_1, dt11_1, dt13_1 
      double precision dt15_1, dt3_2, dt10_2, dt17_2, dt19_2, dt1_3
      double precision dt4_3, dt6_3, dt11_3, dt13_3, dt3_4, dt8_4
      double precision dt17_4, dt1_5, dt4_5, dt11_5
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
      double precision bin11, bin12, bin13, bin14, bin15
      double precision bin16, bin17, bin18, bin19, bin20, bin21
*rak:      integer ls     ! [fixed at 5] angular momentum of block
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=5
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_sh_a: spcart not initialized properly',
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
c  spint(-5) = cint(2)*dtrans(2,-5,5) + cint(7)*dtrans(7,-5,5) + cint(16)*dtrans(16,-5,5)
c  spint(-4) = cint(5)*dtrans(5,-4,5) + cint(12)*dtrans(12,-4,5)
c  spint(-3) = cint(2)*dtrans(2,-3,5) + cint(7)*dtrans(7,-3,5) + cint(9)*dtrans(9,-3,5)
c            + cint(16)*dtrans(16,-3,5) + cint(18)*dtrans(18,-3,5)
c  spint(-2) = cint(5)*dtrans(5,-2,5) + cint(12)*dtrans(12,-2,5) + cint(14)*dtrans(14,-2,5)
c  spint(-1) = cint(2)*dtrans(2,-1,5) + cint(7)*dtrans(7,-1,5) + cint(9)*dtrans(9,-1,5)
c            + cint(16)*dtrans(16,-1,5) + cint(18)*dtrans(18,-1,5) + cint(20)*dtrans(20,-1,5) 
c  spint( 0) = cint(3)*dtrans(3,0,5) + cint(8)*dtrans(8,0,5) + cint(10)*dtrans(10,0,5)
c            + cint(17)*dtrans(17,0,5) + cint(19)*dtrans(19,0,5) + cint(21)*dtrans(21,0,5) 
c  spint( 1) = cint(1)*dtrans(1,1,5)   + cint(4)*dtrans(4,1,5)   + cint(6)*dtrans(6,1,5)
c            + cint(11)*dtrans(11,1,5) + cint(13)*dtrans(13,1,5) + cint(15)*dtrans(15,1,5)
c  spint( 2) = cint(3)*dtrans(3,2,5) + cint(10)*dtrans(10,2,5) + cint(17)*dtrans(17,2,5)
c            + cint(19)*dtrans(19,2,5)
c  spint( 3) = cint(1)*dtrans(1,3,5) + cint(4)*dtrans(4,3,5) + cint(6)*dtrans(6,3,5)
c            + cint(11)*dtrans(11,3,5) + cint(13)*dtrans(13,3,5)
c  spint( 4) = cint(3)*dtrans(3,4,5) + cint(8)*dtrans(8,4,5) + cint(17)*dtrans(17,4,5)
c  spint( 5) = cint(1)*dtrans(1,4,5) + cint(4)*dtrans(4,4,5) + cint(11)*dtrans(11,4,5)
c
      dt2_m5  = dtrans(2,-5,5)
      dt7_m5  = dtrans(7,-5,5)
      dt16_m5 = dtrans(16,-5,5)
      dt5_m4  = dtrans(5,-4,5)     
      dt12_m4 = dtrans(12,-4,5)    
      dt2_m3  = dtrans(2,-3,5)     
      dt7_m3  = dtrans(7,-3,5)     
      dt9_m3  = dtrans(9,-3,5)     
      dt16_m3 = dtrans(16,-3,5)
      dt18_m3 = dtrans(18,-3,5)
      dt5_m2  = dtrans(5,-2,5)     
      dt12_m2 = dtrans(12,-2,5)    
      dt14_m2 = dtrans(14,-2,5)    
      dt2_m1  = dtrans(2,-1,5)
      dt7_m1  = dtrans(7,-1,5)
      dt9_m1  = dtrans(9,-1,5)
      dt16_m1 = dtrans(16,-1,5)
      dt18_m1 = dtrans(18,-1,5)
      dt20_m1 = dtrans(20,-1,5)
      dt3_0   = dtrans(3,0,5)      
      dt8_0   = dtrans(8,0,5)      
      dt10_0  = dtrans(10,0,5)     
      dt17_0  = dtrans(17,0,5)     
      dt19_0  = dtrans(19,0,5)     
      dt21_0  = dtrans(21,0,5)     
      dt1_1   = dtrans(1,1,5)      
      dt4_1   = dtrans(4,1,5)      
      dt6_1   = dtrans(6,1,5)      
      dt11_1  = dtrans(11,1,5)     
      dt13_1  = dtrans(13,1,5)     
      dt15_1  = dtrans(15,1,5)     
      dt3_2   = dtrans(3,2,5)      
      dt10_2  = dtrans(10,2,5)     
      dt17_2  = dtrans(17,2,5)     
      dt19_2  = dtrans(19,2,5)     
      dt1_3   = dtrans(1,3,5)      
      dt4_3   = dtrans(4,3,5)      
      dt6_3   = dtrans(6,3,5)      
      dt11_3  = dtrans(11,3,5)     
      dt13_3  = dtrans(13,3,5)                 
      dt3_4   = dtrans(3,4,5)      
      dt8_4   = dtrans(8,4,5)      
      dt17_4  = dtrans(17,4,5)      
      dt1_5   = dtrans(1,5,5)      
      dt4_5   = dtrans(4,5,5)      
      dt11_5  = dtrans(11,5,5)     
      do i=1,ndima
        do g=1,ngls
          bin1  = blockin( 1,g,i)
          bin2  = blockin( 2,g,i)
          bin3  = blockin( 3,g,i)
          bin4  = blockin( 4,g,i)
          bin5  = blockin( 5,g,i)
          bin6  = blockin( 6,g,i)
          bin7  = blockin( 7,g,i)
          bin8  = blockin( 8,g,i)
          bin9  = blockin( 9,g,i)
          bin10 = blockin(10,g,i)
          bin11 = blockin(11,g,i)
          bin12 = blockin(12,g,i)
          bin13 = blockin(13,g,i)
          bin14 = blockin(14,g,i)
          bin15 = blockin(15,g,i)
          bin16 = blockin(16,g,i)
          bin17 = blockin(17,g,i)
          bin18 = blockin(18,g,i)
          bin19 = blockin(19,g,i)
          bin20 = blockin(20,g,i)
          bin21 = blockin(21,g,i)
c
          blockout(-5,g,i) = bin2*dt2_m5  + bin7*dt7_m5   +
     &                       bin16*dt16_m5
          blockout(-4,g,i) = bin5*dt5_m4  + bin12*dt12_m4 
          blockout(-3,g,i) = bin2*dt2_m3  + bin7*dt7_m3   +
     &                       bin9*dt9_m3  + bin16*dt16_m3 +
     &                       bin18*dt18_m3
          blockout(-2,g,i) = bin5*dt5_m2  + bin12*dt12_m2 +
     &                       bin14*dt14_m2
          blockout(-1,g,i) = bin2*dt2_m1  + bin7*dt7_m1   +
     &                       bin9*dt9_m1  + bin16*dt16_m1 +
     &                       bin18*dt18_m1+ bin20*dt20_m1
          blockout( 0,g,i) = bin3*dt3_0   + bin8*dt8_0    +
     &                       bin10*dt10_0 + bin17*dt17_0  +
     &                       bin19*dt19_0 + bin21*dt21_0
          blockout( 1,g,i) = bin1*dt1_1   + bin4*dt4_1    +
     &                       bin6*dt6_1   + bin11*dt11_1  +
     &                       bin13*dt13_1 + bin15*dt15_1
          blockout( 2,g,i) = bin3*dt3_2   + bin10*dt10_2  +
     &                       bin17*dt17_2 + bin19*dt19_2
          blockout( 3,g,i) = bin1*dt1_3   + bin4*dt4_3    +
     &                       bin6*dt6_3   + bin11*dt11_3  +
     &                       bin13*dt13_3
          blockout( 4,g,i) = bin3*dt3_4   + bin8*dt8_4    +
     &                       bin17*dt17_4
          blockout( 5,g,i) = bin1*dt1_5   + bin4*dt4_5    +
     &                       bin11*dt11_5
        enddo
      enddo
      end
*.......................................................................
      subroutine spcart_a_sh_b(blockin, blockout,
     &    ndima, ndimb, ngls, in_place, print_info)
      implicit none
c
c  transforms a block of integrals where the g function is the middle dimension
c  e.g., blockin(ndima,21,ndimb) and blockout(ndima,11,ndimb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "spcartP.fh"
c
c: passed
      integer ndima  ! [input] leading dimension of block
      integer ndimb  ! [input] trailing dimension of block
      integer ngls   ! [input] general contraction length of ls block
      double precision blockin (ndima,21,ngls,ndimb)   ! [input] matrix 
      double precision blockout(ndima,-5:5,ngls,ndimb) ! [output] matrix 
      logical in_place  ! [input] true if blockin and block out are the same pointer
      logical print_info ! [input] print info
c: local
      integer i,j,g
      double precision dt2_m5, dt7_m5, dt16_m5, dt5_m4, dt12_m4, dt2_m3
      double precision dt7_m3, dt9_m3, dt16_m3, dt18_m3, dt5_m2, dt12_m2
      double precision dt14_m2, dt2_m1, dt7_m1, dt9_m1, dt16_m1, dt18_m1
      double precision dt20_m1, dt3_0, dt8_0, dt10_0, dt17_0, dt19_0
      double precision dt21_0, dt1_1, dt4_1, dt6_1, dt11_1, dt13_1 
      double precision dt15_1, dt3_2, dt10_2, dt17_2, dt19_2, dt1_3
      double precision dt4_3, dt6_3, dt11_3, dt13_3, dt3_4, dt8_4
      double precision dt17_4, dt1_5, dt4_5, dt11_5
      double precision bin1, bin2, bin3, bin4, bin5, bin6 
      double precision bin7, bin8, bin9, bin10
      double precision bin11, bin12, bin13, bin14, bin15
      double precision bin16, bin17, bin18, bin19, bin20, bin21
*rak:      integer ls     ! [fixed at 5] angular momentum of block
c::statement function ----- start
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c::statement function ----- end
c
*rak:      ls=5
*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:      write(6,*)' trans matrix used '
*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
c
      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
     &    ('spcart_a_sh_b: spcart not initialized properly',
     &    sph_cart_init, UNKNOWN_ERR)
c
      if (in_place) then
        write (6,*)' in place transformations are not ready yet '
      endif
c
c  spint(-5) = cint(2)*dtrans(2,-5,5) + cint(7)*dtrans(7,-5,5) + cint(16)*dtrans(16,-5,5)
c  spint(-4) = cint(5)*dtrans(5,-4,5) + cint(12)*dtrans(12,-4,5)
c  spint(-3) = cint(2)*dtrans(2,-3,5) + cint(7)*dtrans(7,-3,5) + cint(9)*dtrans(9,-3,5)
c            + cint(16)*dtrans(16,-3,5) + cint(18)*dtrans(18,-3,5)
c  spint(-2) = cint(5)*dtrans(5,-2,5) + cint(12)*dtrans(12,-2,5) + cint(14)*dtrans(14,-2,5)
c  spint(-1) = cint(2)*dtrans(2,-1,5) + cint(7)*dtrans(7,-1,5) + cint(9)*dtrans(9,-1,5)
c            + cint(16)*dtrans(16,-1,5) + cint(18)*dtrans(18,-1,5) + cint(20)*dtrans(20,-1,5) 
c  spint( 0) = cint(3)*dtrans(3,0,5) + cint(8)*dtrans(8,0,5) + cint(10)*dtrans(10,0,5)
c            + cint(17)*dtrans(17,0,5) + cint(19)*dtrans(19,0,5) + cint(21)*dtrans(21,0,5) 
c  spint( 1) = cint(1)*dtrans(1,1,5)   + cint(4)*dtrans(4,1,5)   + cint(6)*dtrans(6,1,5)
c            + cint(11)*dtrans(11,1,5) + cint(13)*dtrans(13,1,5) + cint(15)*dtrans(15,1,5)
c  spint( 2) = cint(3)*dtrans(3,2,5) + cint(10)*dtrans(10,2,5) + cint(17)*dtrans(17,2,5)
c            + cint(19)*dtrans(19,2,5)
c  spint( 3) = cint(1)*dtrans(1,3,5) + cint(4)*dtrans(4,3,5) + cint(6)*dtrans(6,3,5)
c            + cint(11)*dtrans(11,3,5) + cint(13)*dtrans(13,3,5)
c  spint( 4) = cint(3)*dtrans(3,4,5) + cint(8)*dtrans(8,4,5) + cint(17)*dtrans(17,4,5)
c  spint( 5) = cint(1)*dtrans(1,4,5) + cint(4)*dtrans(4,4,5) + cint(11)*dtrans(11,4,5)
c
      dt2_m5  = dtrans(2,-5,5)     
      dt7_m5  = dtrans(7,-5,5)     
      dt16_m5 = dtrans(16,-5,5)
      dt5_m4  = dtrans(5,-4,5)     
      dt12_m4 = dtrans(12,-4,5)    
      dt2_m3  = dtrans(2,-3,5)     
      dt7_m3  = dtrans(7,-3,5)     
      dt9_m3  = dtrans(9,-3,5)     
      dt16_m3 = dtrans(16,-3,5)
      dt18_m3 = dtrans(18,-3,5)
      dt5_m2  = dtrans(5,-2,5)     
      dt12_m2 = dtrans(12,-2,5)    
      dt14_m2 = dtrans(14,-2,5)    
      dt2_m1  = dtrans(2,-1,5)
      dt7_m1  = dtrans(7,-1,5)
      dt9_m1  = dtrans(9,-1,5)
      dt16_m1 = dtrans(16,-1,5)
      dt18_m1 = dtrans(18,-1,5)
      dt20_m1 = dtrans(20,-1,5)
      dt3_0   = dtrans(3,0,5)      
      dt8_0   = dtrans(8,0,5)      
      dt10_0  = dtrans(10,0,5)     
      dt17_0  = dtrans(17,0,5)     
      dt19_0  = dtrans(19,0,5)     
      dt21_0  = dtrans(21,0,5)     
      dt1_1   = dtrans(1,1,5)      
      dt4_1   = dtrans(4,1,5)      
      dt6_1   = dtrans(6,1,5)      
      dt11_1  = dtrans(11,1,5)     
      dt13_1  = dtrans(13,1,5)     
      dt15_1  = dtrans(15,1,5)     
      dt3_2   = dtrans(3,2,5)      
      dt10_2  = dtrans(10,2,5)     
      dt17_2  = dtrans(17,2,5)     
      dt19_2  = dtrans(19,2,5)     
      dt1_3   = dtrans(1,3,5)      
      dt4_3   = dtrans(4,3,5)      
      dt6_3   = dtrans(6,3,5)      
      dt11_3  = dtrans(11,3,5)     
      dt13_3  = dtrans(13,3,5)                 
      dt3_4   = dtrans(3,4,5)      
      dt8_4   = dtrans(8,4,5)      
      dt17_4  = dtrans(17,4,5)      
      dt1_5   = dtrans(1,5,5)      
      dt4_5   = dtrans(4,5,5)      
      dt11_5  = dtrans(11,5,5)     
      do j=1,ndimb
        do g=1,ngls
          do i=1,ndima
            bin1  = blockin(i, 1,g,j)
            bin2  = blockin(i, 2,g,j)
            bin3  = blockin(i, 3,g,j)
            bin4  = blockin(i, 4,g,j)
            bin5  = blockin(i, 5,g,j)
            bin6  = blockin(i, 6,g,j)
            bin7  = blockin(i, 7,g,j)
            bin8  = blockin(i, 8,g,j)
            bin9  = blockin(i, 9,g,j)
            bin10 = blockin(i,10,g,j)
            bin11 = blockin(i,11,g,j)
            bin12 = blockin(i,12,g,j)
            bin13 = blockin(i,13,g,j)
            bin14 = blockin(i,14,g,j)
            bin15 = blockin(i,15,g,j)
            bin16 = blockin(i,16,g,j)
            bin17 = blockin(i,17,g,j)
            bin18 = blockin(i,18,g,j)
            bin19 = blockin(i,19,g,j)
            bin20 = blockin(i,20,g,j)
            bin21 = blockin(i,21,g,j)
c
            blockout(i,-5,g,j) = bin2*dt2_m5  + bin7*dt7_m5   +
     &                           bin16*dt16_m5
            blockout(i,-4,g,j) = bin5*dt5_m4  + bin12*dt12_m4 
            blockout(i,-3,g,j) = bin2*dt2_m3  + bin7*dt7_m3   +
     &                           bin9*dt9_m3  + bin16*dt16_m3 +
     &                           bin18*dt18_m3
            blockout(i,-2,g,j) = bin5*dt5_m2  + bin12*dt12_m2 +
     &                           bin14*dt14_m2
            blockout(i,-1,g,j) = bin2*dt2_m1  + bin7*dt7_m1   +
     &                           bin9*dt9_m1  + bin16*dt16_m1 +
     &                           bin18*dt18_m1+ bin20*dt20_m1
            blockout(i, 0,g,j) = bin3*dt3_0   + bin8*dt8_0    +
     &                           bin10*dt10_0 + bin17*dt17_0  +
     &                           bin19*dt19_0 + bin21*dt21_0
            blockout(i, 1,g,j) = bin1*dt1_1   + bin4*dt4_1    +
     &                           bin6*dt6_1   + bin11*dt11_1  +
     &                           bin13*dt13_1 + bin15*dt15_1
            blockout(i, 2,g,j) = bin3*dt3_2   + bin10*dt10_2  +
     &                           bin17*dt17_2 + bin19*dt19_2
            blockout(i, 3,g,j) = bin1*dt1_3   + bin4*dt4_3    +
     &                           bin6*dt6_3   + bin11*dt11_3  +
     &                           bin13*dt13_3
            blockout(i, 4,g,j) = bin3*dt3_4   + bin8*dt8_4    +
     &                           bin17*dt17_4
            blockout(i, 5,g,j) = bin1*dt1_5   + bin4*dt4_5    +
     &                           bin11*dt11_5
          enddo
        enddo
      enddo
      end
*.......................................................................
C>
C> \brief Transform 1-electron integrals from Cartesian to spherical
C> harmonic basis functions
C>
C> Even if a basis set is specified as spherical harmonic in the input
C> the code actually calculates the integrals in terms of Cartesian
C> functions. Hence an additional transformation step is needed to 
C> obtain the specified integrals in spherical harmonics. This routine
C> controls that transformation for various angular momenta.
C>
      subroutine spcart_tran1e(
     &    buf, scr,
     &    nbf_xr,nbf_xc,type_r,ngen_r,
     &    nbf_sr,nbf_sc,type_c,ngen_c,
     &    print)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
c
c
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
c
c Can use a special call setting type_? to 0 for basis sets that 
c are mixed   i_bas is spherical and j_bas is not dimensions 
c are not based on type only actions.
c
c      
      integer nbf_xr !< [input] number of rows of cartesian block
      integer nbf_xc !< [input] number of columns of cartesian block
      integer nbf_sr !< [input] number of rows of spherical block
      integer nbf_sc !< [input] number of columns of spherical block
      integer ngen_r !< [input] general contraction length for rows
      integer ngen_c !< [input] general contraction length for columns
      integer type_r !< [input] angular momentem for rows
      integer type_c !< [input] angular momentem for columns
      double precision buf(*) !< [input/output] cartesian block on input
                              !< and spherical block on output      
      double precision scr(*) !< [scratch] use to hold half transformed block
      logical print           !< [input] print integrals at each stage of the
                              !< transformation (cart/half/spherical)
c:: local
      integer ngen_rc
      logical problem_sp
      ngen_rc = ngen_r*ngen_c
c... more error checking
      problem_sp = (type_r.eq.-1).and.(ngen_r.gt.1)
      problem_sp = problem_sp.or.((type_c.eq.-1).and.(ngen_c.gt.1))
      if (problem_sp) then
        write(6,*)' spcart_tran1e: sp function error '
        write(6,*)' sp function block passed with more ',
     &      'than one general contraction: ',
     &      'this is not allowed in NWChem'
        write(6,*)' type r',type_r
        write(6,*)' ngen r',ngen_r
        write(6,*)' type c',type_c
        write(6,*)' ngen c',ngen_c
        call errquit('spcart_tran1e: fatal error',911, UNKNOWN_ERR)
      endif
c
      if (type_r.lt.2.and.type_c.lt.2) then
c.................................. neither c or r need to be transformed 
c                                   (X,Y) X is s, p, or l and Y is s, p, or l
        if (print) then
          write(luout,*)
     &        ' cartesian matrix and spherical matrix the same '
          call output(buf,1,nbf_xr*ngen_r,1,nbf_xc*ngen_c,
     &          nbf_xr*ngen_r,nbf_xc*ngen_c,1)
        endif
c
      elseif (type_r.lt.2.and.type_c.ge.2) then
c.................................. c needs to be transformed
* print cartesian matrix  
* buf is buf(spherical,cartesian)
        if (print) then
          write(luout,*)' cartesian matrix '
          call output(buf,1,nbf_xr*ngen_r,1,nbf_xc*ngen_c,
     &          nbf_xr*ngen_r,nbf_xc*ngen_c,1)
        endif
* scr is buf(spherical,cartesian) ! copy it
        call dcopy((nbf_xr*nbf_xc*ngen_rc),buf,1,scr,1)
        if (nbf_xr.ne.nbf_sr) call errquit
     &      ('spcart_tran1e: nbf_xr.ne.nbf_sr  (xr-sr) =',
     &      (nbf_xr-nbf_sr), UNKNOWN_ERR)
        call spcart_a_s(scr,buf,(ngen_r*nbf_sr),type_c,ngen_c,
     &      .false.,print)
        if (print) then
          write(luout,*)' spherical matrix '
          call output(buf,1,ngen_r*nbf_sr,1,ngen_c*nbf_sc,
     &          ngen_r*nbf_sr,ngen_c*nbf_sc,1)
        endif
      elseif (type_r.ge.2.and.type_c.lt.2) then
c.................................. r needs to be transformed
* print cartesian matrix  
* buf is buf(cartesian,cartesian)
        if (print) then
          write(luout,*)' cartesian matrix '
          call output(buf,1,ngen_r*nbf_xr,1,ngen_c*nbf_xc,
     &          ngen_r*nbf_xr,ngen_c*nbf_xc,1)
        endif
* scr is buf(spherical,cartesian) ! copy it
        call dcopy((ngen_rc*nbf_xr*nbf_xc),buf,1,scr,1)
        if (nbf_xc.ne.nbf_sc) call errquit
     &      ('spcart_tran1e: nbf_xc.ne.nbf_sc  (xc-sc) =',
     &      (nbf_xc-nbf_sc), UNKNOWN_ERR)
        call spcart_s_a(scr,buf,ngen_c*nbf_sc,type_r,ngen_r,
     &      .false.,print)
        if (print) then
          write(luout,*)' spherical matrix '
          call output(buf,1,ngen_r*nbf_sr,1,ngen_c*nbf_sc,
     &          ngen_r*nbf_sr,ngen_c*nbf_sc,1)
        endif
      elseif (type_r.ge.2.and.type_c.ge.2) then
c.................................. both r and c need to be transformed
* print cartesian matrix  
* buf is buf(cartesian,cartesian)
        if (print) then
          write(luout,*)' cartesian matrix '
          call output(buf,1,ngen_r*nbf_xr,1,ngen_c*nbf_xc,
     &          ngen_r*nbf_xr,ngen_c*nbf_xc,1)
        endif
*
*... buf(xr,xc) -> scr(xr,sc) : scr is half transformed matrix
        call spcart_a_s(buf,scr,ngen_r*nbf_xr,type_c,ngen_c,
     &      .false.,print)
* print half transformed matrix
        if (print) then
          write(luout,*)' half cartesian half spherical matrix '
          call output(scr,1,ngen_r*nbf_xr,1,ngen_c*nbf_sc,
     &          ngen_r*nbf_xr,ngen_c*nbf_sc,1)
        endif
*
*... scr(xr,sc) -> buf(sr,sc)
        call spcart_s_a(scr,buf,ngen_c*nbf_sc,type_r,ngen_r,
     &      .false.,print)
* print spherical block
        if (print) then
          write(luout,*)' spherical matrix '
          call output(buf,1,ngen_r*nbf_sr,1,ngen_c*nbf_sc,
     &          ngen_r*nbf_sr,ngen_c*nbf_sc,1)
        endif
      else
        write(luout,*) ' case not possible '
        call errquit('spcart_tran1e: should never get here ? ',911,
     &       UNKNOWN_ERR)
      endif
c
      end
*.......................................................................
      subroutine spcart_bra2etran(
     &    buf, scr,
     &    nbf_xj,nbf_xi,
     &    nbf_sj,nbf_si,
     &    type_j,type_i,
     &    ngen_j,ngen_i,
     &    ndim_ket,
     &    print)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
c
c routine that transforms a 2e cartesian block buf_cart(ndim_ket,nbf_xj,nbf_xi) to 
c    a spherical block buf_sph(ndim_ket,nbf_sj,nbf_si) 
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
      integer ngen_j, ngen_i  ! [input] general contraction length for j and i
      integer ndim_ket
      double precision buf(*) ! [input/output] cartesian block on input
*.............................!   and spherical block on output      
      double precision scr(*) ! [scratch] use to hold half transformed block
      logical print           ! [input} print integrals at each stage of the
*.............................!   transformation (cart/half/spherical)
*::local
      logical problem_sp
      integer ngen_ij, sxi, sxj, ssi, ssj
c
      if(print) then
         write(6,*)' bra2etran '
         write(6,*)'bra2: nbf_xj, nbf_xi :',nbf_xj, nbf_xi
         write(6,*)'bra2: nbf_sj, nbf_si :',nbf_sj, nbf_si
         write(6,*)'bra2: type_j, type_i :',type_j, type_i
         write(6,*)'bra2: ngen_j, ngen_i :',ngen_j, ngen_i
         write(6,*)'bra2: ndim_ket       :',ndim_ket
      endif
c... more error checking
      problem_sp = (type_j.eq.-1).and.(ngen_j.gt.1)
      problem_sp = problem_sp.or.((type_i.eq.-1).and.(ngen_i.gt.1))
      if (problem_sp) then
        write(6,*)' spcart_bra2etran: sp function error '
        write(6,*)' sp function block passed with more ',
     &      'than one general contraction: ',
     &      'this is not allowed in NWChem'
        write(6,*)' type j',type_j
        write(6,*)' ngen j',ngen_j
        write(6,*)' type i',type_i
        write(6,*)' ngen i',ngen_i
        call errquit('spcart_bra2etran: fatal error',911, UNKNOWN_ERR)
      endif
c
      ngen_ij = ngen_i*ngen_j
      sxi = nbf_xi*ngen_i
      sxj = nbf_xj*ngen_j
      ssi = nbf_si*ngen_i
      ssj = nbf_sj*ngen_j
c
      if (type_j.lt.2.and.type_i.lt.2) then
c.................................. neither i or j need to be transformed 
c                                   (X,Y) X is s, p, or l and Y is s, p, or l
        if (print) then
          write(luout,*)
     &        ' cartesian matrix and spherical matrix the same '
          call slice_am_print(ndim_ket,sxj,sxi,buf,
     &        ' (ket:j-cart/s:i-cart/s')
        endif
      elseif (type_j.lt.2.and.type_i.ge.2) then
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c.................................. i needs to be transformed
* print cartesian matrix  
* buf is buf(ndim_ket,j_spherical,i_cartesian)
        if (print) then
          write(luout,*)' (ket:j-spherical:i-cartesian) matrix '
          call slice_am_print(ndim_ket,ssj,sxi,buf,
     &          ' (ket:j-spherical:i-cartesian) matrix ')
        endif
* scr is buf(ndim_ket,j_spherical,i_cartesian) ! copy it
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
        call dcopy((ndim_ket*ssj*sxi),buf,1,scr,1)
        if (nbf_xj.ne.nbf_sj) call errquit
     &      ('spcart_bra2etran: nbf_xj.ne.nbf_sj  (xj-sj) =',
     &      (nbf_xj-nbf_sj), UNKNOWN_ERR)
        call spcart_a_s(scr,buf,(ndim_ket*ssj),
     &        type_i,ngen_i,.false.,print)
        if (print) then
          write(luout,*)' (ket:j-spherical:i-shperical) matrix '
          call slice_am_print(ndim_ket,ssj,ssi,buf,
     &          ' (ket:j-spherical:i-shperical) matrix ')
        endif
      elseif (type_j.ge.2.and.type_i.lt.2) then
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c.................................. j needs to be transformed
* print cartesian matrix  
* buf is buf(ndim_ket,j_cartesian,i_spherical)
        if (print) then
          write(luout,*)' (ket:j-cartesian:i-spherical) matrix '
          call slice_am_print(ndim_ket,sxj,ssi,buf,
     &          ' (ket:j-cartesian:i-spherical) matrix ')
        endif
* buf is buf(ndim_ket,j_cartesian,i_spherical)
* scr is buf(ndim_ket,j_cartesian,i_spherical) ! copy it
        call dcopy((ndim_ket*sxj*ssi),buf,1,scr,1)
        if (nbf_xi.ne.nbf_si) call errquit
     &      ('spcart_bra2etran: nbf_xc.ne.nbf_sc  (xi-si) =',
     &      (nbf_xi-nbf_si), UNKNOWN_ERR)
        call spcart_a_s_b(scr,buf,ndim_ket,ssi,
     &        type_j,ngen_j,.false.,print)
        if (print) then
          write(luout,*)' (ket:j-spherical:i-spherical) matrix '
          call slice_am_print(ndim_ket,ssj,ssi,buf,
     &          ' (ket:j-spherical:i-spherical) matrix ')
        endif
      elseif (type_j.ge.2.and.type_i.ge.2) then
*...................ERI(ndim_ket,jbf_lo:jbf_hi,ibf_lo,ibf_hi)
c........................ both j and i need to be transformed
* print cartesian matrix  
* buf is buf(ndim_ket,j_cartesian,i_cartesian)
        if (print) then
          write(luout,*)' (ket:j-cartesian:i-cartesian) matrix '
          call slice_am_print(ndim_ket,sxj,sxi,buf,
     &          ' (ket:j-cartesian:i-cartesian) matrix ')
        endif
*
*... buf(ndim_ket,xj,xi) -> scr(xj,si) : scr is half transformed matrix
        call spcart_a_s(buf,scr,(ndim_ket*sxj),
     &        type_i,ngen_i,.false.,print)
* print half transformed matrix
        if (print) then
          write(luout,*)' (ket:j-cartesian:i-spherical) matrix '
          call slice_am_print(ndim_ket,sxj,ssi,scr,
     &          ' (ket:j-cartesian:i-spherical) matrix ')
        endif
*
*... scr(xj,si) -> buf(sj,si)
        call spcart_a_s_b(scr,buf,ndim_ket,ssi,
     &        type_j,ngen_j,.false.,print)
* print spherical block
        if (print) then
          write(luout,*)' (ket:j-spherical:i-spherical) matrix '
          call slice_am_print(ndim_ket,ssj,ssi,buf,
     &          ' (ket:j-spherical:i-spherical) matrix ')
        endif
c
      else
        write(luout,*) ' case not possible '
        call errquit('spcart_bra2etran: should never get here ? ',911,
     &       UNKNOWN_ERR)
      endif
c
      end
*.......................................................................
      subroutine spcart_ket2etran(
     &    buf, scr,
     &    nbf_xl,nbf_xk,
     &    nbf_sl,nbf_sk,
     &    type_l,type_k,
     &    ngen_l,ngen_k,
     &    ndim_bra,
     &    print)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
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
      integer type_l, type_k  ! [input] angular momentem for l and k
      integer ngen_l, ngen_k  ! [input] general contraction length for l and k
      integer ndim_bra
      double precision buf(*) ! [input/output] cartesian block on input
*.............................!   and spherical block on output      
      double precision scr(*) ! [scratch] use to hold half transformed block
      logical print           ! [input} print integrals at each stage of the
*.............................!   transformation (cart/half/spherical)
*::local
      logical problem_sp
      integer ngen_kl, sxl, sxk, ssl, ssk
c
      if(print) then
         write(6,*)' ket2etran '
         write(6,*)'ket2: nbf_xl, nbf_xk :',nbf_xl, nbf_xk
         write(6,*)'ket2: nbf_sl, nbf_sk :',nbf_sl, nbf_sk
         write(6,*)'ket2: type_l, type_k :',type_l, type_k
         write(6,*)'ket2: ngen_l, ngen_k :',ngen_l, ngen_k
         write(6,*)'ket2: ndim_bra       :',ndim_bra
      endif
c... more error checking
      problem_sp = (type_l.eq.-1).and.(ngen_l.gt.1)
      problem_sp = problem_sp.or.((type_k.eq.-1).and.(ngen_k.gt.1))
      if (problem_sp) then
        write(6,*)' spcart_ket2etran: sp function error '
        write(6,*)' sp function block passed with more ',
     &      'than one general contraction: ',
     &      'this is not allowed in NWChem'
        write(6,*)' type l',type_l
        write(6,*)' ngen l',ngen_l
        write(6,*)' type k',type_k
        write(6,*)' ngen k',ngen_k
        call errquit('spcart_ket2etran: fatal error',911,
     &       UNKNOWN_ERR)
      endif
c
c
      ngen_kl = ngen_l * ngen_k
      sxl = nbf_xl*ngen_l
      sxk = nbf_xk*ngen_k
      ssl = nbf_sl*ngen_l
      ssk = nbf_sk*ngen_k
c
      if (type_l.lt.2.and.type_k.lt.2) then
c.................................. neither k or l need to be transformed 
c                                   (X,Y) X is s, p, or l and Y is s, p, or l
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
        if (print) then
          write(luout,*)
     &        ' cartesian matrix and spherical matrix the same '
          call slice_ma_print(ndim_bra,sxl,sxk,buf,
     &          ' cartesian matrix and spherical matrix the same ')
        endif
      elseif (type_l.lt.2.and.type_k.ge.2) then
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
c.................................. k needs to be transformed
* print cartesian matrix  
* buf is buf(spherical,cartesian,ndim_bra)
        if (print) then
          write(luout,*)' (spherical:cartesian:bra) matrix '
          call slice_ma_print(ndim_bra,ssl,sxk,buf,
     &          ' (l-spherical:k-cartesian:bra) matrix ')
        endif
* scr is buf(l_spherical,k_cartesian,ndim_bra) ! copy it
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
        call dcopy((ndim_bra*ssl*sxk),buf,1,scr,1)
        if (nbf_xl.ne.nbf_sl) call errquit
     &      ('spcart_ket2etran: nbf_xl.ne.nbf_sl  (xl-sl) =',
     &      (nbf_xl-nbf_sl), UNKNOWN_ERR)
        call spcart_a_s_b(scr,buf,ssl,ndim_bra,
     &        type_k,ngen_k,.false.,print)
        if (print) then
          write(luout,*)' (l-spherical:k-shperical:bra) matrix '
          call slice_ma_print(ndim_bra,ssl,ssk,buf,
     &          ' (l-spherical:k-shperical:bra) matrix ')
        endif
      elseif (type_l.ge.2.and.type_k.lt.2) then
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
c.................................. l needs to be transformed
* print cartesian matrix  
* buf is buf(l_cartesian,k_spherical,ndim_bra)
        if (print) then
          write(luout,*)' (l-cartesian:k-spherical:bra) matrix '
          call slice_ma_print(ndim_bra,sxl,ssk,buf,
     &          ' (l-cartesian:k-spherical:bra) matrix ')
        endif
* buf is buf(l_cartesian,k_spherical,ndim_bra)
* scr is buf(l_cartesian,k_spherical,ndim_bra)
        call dcopy((ndim_bra*sxl*ssk),buf,1,scr,1)
        if (nbf_xk.ne.nbf_sk) call errquit
     &      ('spcart_ket2etran: nbf_xk.ne.nbf_sk  (xk-sk) =',
     &      (nbf_xk-nbf_sk), UNKNOWN_ERR)
        call spcart_s_a(scr,buf,(ndim_bra*ssk),
     &        type_l,ngen_l,.false.,print)
        if (print) then
          write(luout,*)' (l-spherical:k-spherical:bra) matrix '
          call slice_ma_print(ndim_bra,ssl,ssk,buf,
     &          ' (l-spherical:k-spherical:bra) matrix ')
        endif
      elseif (type_l.ge.2.and.type_k.ge.2) then
*...................ERI(lbf_lo:lbf_hi,kbf_lo,kbf_hi,ndim_bra)
c........................ both k and l need to be transformed
* print cartesian matrix  
* buf is buf(l_cartesian,k_cartesian,ndim_bra)
        if (print) then
          write(luout,*)' (l-cartesian:k-cartesian:bra) matrix '
          call slice_ma_print(ndim_bra,sxl,sxk,buf,
     &          ' (l-cartesian:k-cartesian:bra) matrix ')
        endif
*
*... buf(xl,xk) -> scr(sl,xk,ndim_bra,) : scr is half transformed matrix
        call spcart_s_a(buf,scr,(sxk*ndim_bra),
     &        type_l,ngen_l,.false.,print)
* print half transformed matrix
        if (print) then
          write(luout,*)' (l-spherical:k-cartesian:bra) matrix '
          call slice_ma_print(ndim_bra,sxl,ssk,scr,
     &          ' (l-cartesian:k-spherical:bra) matrix ')
        endif
*
*... scr(sl,xk) -> buf(sl,sk)
        call spcart_a_s_b(scr,buf,ssl,ndim_bra,
     &        type_k,ngen_k,.false.,print)
* print spherical block
        if (print) then
          write(luout,*)' (l-spherical:k-spherical:bra) matrix '
          call slice_ma_print(ndim_bra,ssl,ssk,buf,
     &          ' (l-spherical:k-spherical:bra) matrix ')
        endif
c
      else
        write(luout,*) ' case not possible '
        call errquit('spcart_ket2etran: should never get here ? ',911,
     &       UNKNOWN_ERR)
      endif
c
      end
      subroutine slice_am_print(na,ni,nj,M,msg)
      implicit none
*
* routine to print 2d slices(i,j) of a matrix dimensioned
*
*  Matrix(na,ni,nj)
*
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"      
      integer na, ni, nj           ! [input] matrix dimensions
      double precision M(na,ni,nj) ! [input] matrix to print
      character*(*) msg            ! [input] message for output
c
      integer h_slice, k_slice 
      integer a,i,j,cnt
c
      if (.not.ma_push_get(mt_dbl,(ni*nj),
     &      'slice_am_print scratch',h_slice,k_slice))
     &      call errquit
     &      ('slice_am_print: could not allocate (push) slice',911,
     &      UNKNOWN_ERR)
c
      write(luout,*)' sliced matrix output ',msg
      write(luout,*)na,' slices to be printed ',msg
      do a=1,na
        call dcopy((ni*nj),0.0d00,0,dbl_mb(k_slice),1)
        cnt = k_slice
        do i=1,ni
          do j=1,nj
            dbl_mb(cnt) = M(a,i,j)
            cnt = cnt + 1
          enddo
        enddo
        write(luout,*)' slice ',a,' ',msg
        call output(dbl_mb(k_slice),1,ni,1,nj,ni,nj,1)
      enddo
c
      if (.not.ma_pop_stack(h_slice))call errquit
     &      ('slice_am_print: could not deallocate (pop) slice',911,
     &      MEM_ERR)
c
      end
      subroutine slice_ma_print(na,ni,nj,M,msg)
      implicit none
*
* routine to print 2d slices(i,j) of a matrix dimensioned
*
*  Matrix(ni,nj,na)
*
#include "stdio.fh"
      integer na, ni, nj           ! [input] matrix dimensions
      double precision M(ni,nj,na) ! [input] matrix to print
      character*(*) msg            ! [input] message for output
c
      integer a
      write(luout,*)' sliced matrix output ',msg
      write(luout,*)na,' slices to be printed ',msg
      do a=1,na
        write(luout,*)' slice ',a,' ',msg
        call output(M(1,1,a),1,ni,1,nj,ni,nj,1)
      enddo
      end
      subroutine slice_amb_print(na,nb,ni,nj,M,msg)
*
* routine to print 2d slices(i,j) of a matrix dimensioned
*
*  Matrix(na,ni,nj,nb)
*
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"      
      integer na, nb, ni, nj          ! [input] matrix dimensions
      double precision M(na,ni,nj,nb) ! [input] matrix to print
      character*(*) msg               ! [input] message for output
c
      integer h_slice, k_slice 
      integer a,i,j,b,cnt
c
      if (.not.ma_push_get(mt_dbl,(ni*nj),
     &      'slice_amb_print scratch',h_slice,k_slice))
     &      call errquit
     &      ('slice_amb_print: could not allocate (push) slice',911,
     &      MEM_ERR)
c
      write(luout,*)' sliced matrix output ',msg
      write(luout,*)(na*nb),' slices to be printed ',msg
      do a=1,na
        do b=1,nb
          call dcopy((ni*nj),0.0d00,0,dbl_mb(k_slice),1)
          cnt = k_slice
          do i=1,ni
            do j=1,nj
              dbl_mb(cnt) = M(a,i,j,b)
              cnt = cnt + 1
            enddo
          enddo
          write(luout,*)' slice ',a,b,' ',msg
          call output(dbl_mb(k_slice),1,ni,1,nj,ni,nj,1)
        enddo
      enddo
c
      if (.not.ma_pop_stack(h_slice))call errquit
     &      ('slice_amb_print: could not deallocate (pop) slice',911,
     &      MEM_ERR)
c
      end
      subroutine spcart_cart_overlap(lval,rankov,overlap)
      implicit none
*
* compute the cartsian overlap matrix for a given type of basis function
* This matrix is used to form the inverse of dtrans
*
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "basdeclsP.fh"
      integer lval
      integer rankov
      double precision overlap(rankov,rankov)
*      
      double precision double_dummy
      double precision xyz(3)
      double precision zeta, coef
      integer h_scr,k_scr
      integer scr_size
*
*      write(6,*)' spcart_cart_overlap verify 1',' lval   = ',lval,
*     &      ' rankov = ',rankov
*      if (.not. ma_verify_allocator_stuff())
*     &    stop ' spcart_cart_overlap '
      call dcopy (3,0.0d00,0,xyz,1)  ! do it at the origin
      zeta  = 1.2345d00
      coef = 1.0d00
*
* normalize exponents just like in integrals       
      call nmcoef(zeta,coef,lval,1,BasNorm_STD)
*      
      scr_size = 200000
*      write(6,*)' spcart_cart_overlap verify 2',' lval   = ',lval,
*     &      ' rankov = ',rankov
*      if (.not. ma_verify_allocator_stuff())
*     &    stop ' spcart_cart_overlap '
      call hf1(xyz,zeta,coef,1,1,lval,
     &    xyz,zeta,coef,1,1,lval,
     &    xyz,double_dummy,double_dummy,1,
     &    overlap,double_dummy,double_dummy,
     &    (rankov*rankov),
     &    .true.,.false.,.false.,.false.,
     &    .true.,
     &    double_dummy,scr_size)
      scr_size = max(scr_size,5000)
*      write(6,*)' spcart_cart_overlap verify 3',' lval   = ',lval,
*     &      ' rankov = ',rankov
*      if (.not. ma_verify_allocator_stuff())
*     &    stop ' spcart_cart_overlap '
*      write(6,*)' scr_size is ',scr_size
      if (.not.ma_push_get(mt_dbl,scr_size,
     &    'spcart overlap scratch buffer',
     &    h_scr,k_scr)) call errquit
     &    ('spcart_cart_overlap: could not allocate scratch buffer',
     &    911, MEM_ERR)
      call dcopy(scr_size,0.0d00,0,dbl_mb(k_scr),1)
*      write(6,*)' spcart_cart_overlap verify 4',' lval   = ',lval,
*     &      ' rankov = ',rankov
*      if (.not. ma_verify_allocator_stuff())
*     &    stop ' spcart_cart_overlap '
      call dcopy((rankov*rankov),0.0d00,0,overlap,1)
*      write(6,*)' spcart_cart_overlap verify 5',' lval   = ',lval,
*     &      ' rankov = ',rankov
*      if (.not. ma_verify_allocator_stuff())
*     &    stop ' spcart_cart_overlap '
*      write(6,*)' calling hf1 '
      call hf1(xyz,zeta,coef,1,1,lval,
     &    xyz,zeta,coef,1,1,lval,
     &    xyz,double_dummy,double_dummy,1,
     &    overlap,double_dummy,double_dummy,
     &    (rankov*rankov),
     &    .true.,.false.,.false.,.false.,
     &    .false.,
     &    dbl_mb(k_scr),scr_size)
*      write(6,*)' spcart_cart_overlap verify 6',' lval   = ',lval,
*     &      ' rankov = ',rankov
*      if (.not. ma_verify_allocator_stuff())
*     &    stop ' spcart_cart_overlap '
      if (.not.ma_pop_stack(h_scr))call errquit
     &    ('spcart_cart_overlap: could not pop_stack scratch buffer',
     &    911, MEM_ERR)
*-debug-s
*      write(luout,*)' overlap matrix for lval = ',lval
*      call output(overlap,1,rankov,1,rankov,rankov,rankov,1)
*-debug-s
      end

*.......................................................................
*rak:      subroutine spcart_all(blockin, blockout,li,lj,lk,ll,ni,nj,nk,nl)
*rak:      implicit none
*rak:c
*rak:c  transforms a block of integrals with loop level crap 
*rak:c
*rak:#include "mafdecls.fh"
*rak:#include "errquit.fh"
*rak:#include "spcartP.fh"
*rak:c: functions
*rak:c: passed
*rak:      integer li, lj, lk, ll
*rak:      integer di, dj, dk, dl
*rak:      integer L2i, L2j, L2k, L2l
*rak:      integer ni, nj, nk, nl
*rak:      double precision blockin(
*rak:     &    ((ll+1)*(ll+2)/2)*nl,
*rak:     &    ((lk+1)*(lk+2)/2)*nk,
*rak:     &    ((lj+1)*(lj+2)/2)*nj,
*rak:     &    ((li+1)*(li+2)/2)*ni)
*rak:      double precision blockout(
*rak:     &    (2*ll+1)*nl,
*rak:     &    (2*lk+1)*nk,
*rak:     &    (2*lj+1)*nj,
*rak:     &    (2*li+1)*ni)
*rak:*rak:      double precision blockin(
*rak:*rak:     &    ((ll+1)*(ll+2)/2),nl,
*rak:*rak:     &    ((lk+1)*(lk+2)/2),nk,
*rak:*rak:     &    ((lj+1)*(lj+2)/2),nj,
*rak:*rak:     &    ((li+1)*(li+2)/2),ni)
*rak:*rak:      double precision blockout(
*rak:*rak:     &    -ll:ll,nl,
*rak:*rak:     &    -lk:lk,nk,
*rak:*rak:     &    -lj:lj,nj,
*rak:*rak:     &    -li:li,ni)
*rak:      integer i,j,k,l,ig,jg,kg,lg,is,js,ks,ls
*rak:      integer in_x_i,in_x_j,in_x_k,in_x_l
*rak:      integer in_s_i,in_s_j,in_s_k,in_s_l
*rak:      integer sizeblocks
*rak:c
*rak:      double precision sum, sumadd
*rak:c
*rak:c::statement function ----- start
*rak:      integer iic,iis,iil
*rak:      integer sfi,sfj,sfll
*rak:      integer sf_indx, sf_inds
*rak:      double precision Dtrans
*rak:      Dtrans(iic,iis,iil) =
*rak:     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
*rak:     &           ((iis+iil)*(iil+1)*(iil+2)/2)
*rak:     &           + iic - 1)
*rak:      sf_indx(sfi,sfj,sfll) = (sfj-1)*(sfll+1)*(sfll+2)/2 + sfi
*rak:      sf_inds(sfi,sfj,sfll) = (sfj-1)*(2*sfll+1) + sfi + sfll + 1
*rak:c::statement function ----- end
*rak:*rak:      ls=4
*rak:*rak:      write(6,*)'a_s,(ndima,l2s) ndima = ',ndima,'   ls = ',ls
*rak:*rak:      write(6,*)' trans matrix used '
*rak:*rak:      call output(dbl_mb((int_mb(k_sp2c_lindx+ls))),1,
*rak:*rak:     &    ((ls+1)*(ls+2)/2),1,(2*ls+1),
*rak:*rak:     &    ((ls+1)*(ls+2)/2),(2*ls+1),1)
*rak:c
*rak:      if (sph_cart_init.ne.SPH_CART_INIT_VALUE) call errquit
*rak:     &    ('spcart_a_sg_b: spcart not initialized properly',
*rak:     &    sph_cart_init)
*rak:c
*rak:c
*rak:      write(6,*)' *** Li *** ',Li
*rak:      call spcart_print_dtrans(Li)
*rak:      write(6,*)' *** Lj *** ',Lj
*rak:      call spcart_print_dtrans(Lj)
*rak:      write(6,*)' *** Lk *** ',Lk
*rak:      call spcart_print_dtrans(Lk)
*rak:      write(6,*)' *** Ll *** ',Ll
*rak:      call spcart_print_dtrans(Ll)
*rak:      sizeblocks = ni*(2*li+1)
*rak:      sizeblocks = sizeblocks*nj*(2*lj+1)
*rak:      sizeblocks = sizeblocks*nk*(2*lk+1)
*rak:      sizeblocks = sizeblocks*nl*(2*ll+1)
*rak:      call dfill(sizeblocks,0.0d00,blockout,1)
*rak:      l2i = (li+1)*(li+2)/2
*rak:      l2j = (lj+1)*(lj+2)/2
*rak:      l2k = (lk+1)*(lk+2)/2
*rak:      l2l = (lj+1)*(lj+2)/2
*rak:      do ig = 1,ni
*rak:        do jg = 1,nj
*rak:          do kg = 1,nk
*rak:            do lg = 1,nk
*rak:              do is = -li,li
*rak:                do js = -lj,lj
*rak:                  do ks = -lk,lk
*rak:                    do ls = -ll,ll
*rak:                      sum = 0.0d00
*rak:                      do i = 1,l2i
*rak:                        di = Dtrans(i,is,Li)
*rak:                        do j = 1,l2j
*rak:                          dj = Dtrans(j,js,Lj)
*rak:                          do k = 1,l2k
*rak:                            dk = Dtrans(k,ks,Lk)
*rak:                            do l = 1,l2l
*rak:                              dl = Dtrans(l,ls,Ll)
*rak:*                              sumadd = di*dj*dk*dl*
*rak:*     &                            blockin(l,lg,k,kg,j,jg,i,ig)
*rak:                              in_x_l = sf_indx(l,lg,Ll)
*rak:                              in_x_k = sf_indx(k,kg,Lk)
*rak:                              in_x_j = sf_indx(j,jg,Lj)
*rak:                              in_x_i = sf_indx(i,ig,Li)
*rak:                              sumadd = di*dj*dk*dl*
*rak:     &                            blockin(in_x_l,in_x_k,in_x_j,in_x_i)
*rak:                              sum = sum + sumadd
*rak:                            enddo
*rak:                          enddo
*rak:                        enddo
*rak:                      enddo
*rak:*                      blockout(ls,lg,ks,kg,js,jg,is,ig) = sum
*rak:                      in_s_l = sf_inds(ls,lg,Ll)
*rak:                      in_s_k = sf_inds(ks,kg,Lk)
*rak:                      in_s_j = sf_inds(js,jg,Lj)
*rak:                      in_s_i = sf_inds(is,ig,Li)
*rak:                      blockout(in_s_l,in_s_k,in_s_j,in_s_i) = sum
*rak:                    enddo
*rak:                  enddo
*rak:                enddo
*rak:              enddo
*rak:            enddo
*rak:          enddo
*rak:        enddo
*rak:      enddo
*rak:      end
*.......................................................................
      subroutine spcart_2ctran(buf,scr,lscr,
     &    nbfxi,nbfsi,typei,ngeni,trani,
     &    nbfxj,nbfsj,typej,ngenj,tranj,
     &    printit)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
*
* routine to transform a 1e or 2e 2 center block of integrals
* block(jlo:jhi,:ilo:ihi) 
* 
      double precision buf(*)! [input/output] cartesian/shperical ints.
      integer lscr           ! [input] length of scratch array
      double precision
     &        scr(lscr)      ! [scratch] scratch array
      integer nbfxi          ! [input] no. of cartesian basis funcs:ish
      integer nbfsi          ! [input] no. of spherical basis funcs:ish
      integer typei          ! [input] angular momentum type:ish
      integer ngeni          ! [input] number general contractions:ish
      integer nbfxj          ! [input] no. of cartesian basis funcs:jsh
      integer nbfsj          ! [input] no. of spherical basis funcs:jsh
      integer typej          ! [input] angular momentum type:jsh
      integer ngenj          ! [input] number general contractions:jsh
      logical trani          ! [input] true -> transform ish block
      logical tranj          ! [input] true -> transform jsh block
      logical printit        ! [input] true -> print transform steps
*::local
      integer dimij, dimXX, dimSS
      integer dimj, dimi
      logical problem_sp
      logical FF, FT
      FF = .false.
      FT = .true.
      problem_sp =
     &    (typei.eq.-1).and.(ngeni.gt.1)
      problem_sp = problem_sp.and.
     &    (typej.eq.-1).and.(ngenj.gt.1)
      if (problem_sp) then
        write(luout,*)' spcart_2ctran: sp function error '
        write(luout,*)' sp function block passed with more ',
     &      'than one genereal contraction: ',
     &      'this is not allowed in NWChem'
        write(luout,*)' type i ',typei
        write(luout,*)' ngen i ',ngeni
        write(luout,*)' type j ',typej
        write(luout,*)' ngen j ',ngenj
        call errquit('spcart_2ctran: fatal error',911, UNKNOWN_ERR)
      endif
*
* sanity checking for spherical transforms
*
      if (trani.and.typei.ge.2.and.nbfxi.le.nbfsi) then
        write(luout,*)' sanity check error on i shell info'
        write(luout,*)' shell type     : ',typei
        write(luout,*)' cartesian size :',nbfxi
        write(luout,*)' spherical size :',nbfsi
        call errquit('spcart_2ctran:i: fatal error',911, UNKNOWN_ERR)
      endif
      if (tranj.and.typej.ge.2.and.nbfxj.le.nbfsj) then
        write(luout,*)' sanity check error on j shell info'
        write(luout,*)' shell type     : ',typej
        write(luout,*)' cartesian size :',nbfxj
        write(luout,*)' spherical size :',nbfsj
        call errquit('spcart_2ctran:j: fatal error',911, UNKNOWN_ERR)
      endif
*
* Check for all s, p, sp shells (e.g., no transform required).
*
      dimXX = nbfxi * nbfxj
      dimSS = nbfsi * nbfsj
      if (dimXX.eq.dimSS) return
*
* check scratch size
*
      dimij = nbfxi * ngeni * nbfxj * ngenj
      if (dimij.gt.lscr) then
        write(luout,*)' dimij   :',dimij
        write(luout,*)' lscr    :',lscr
        call errquit
     &      ('spcart_2ctran: error in scratch size for',911, MEM_ERR)
      endif
*
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (trani.and.tranj) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartj,carti) --> scr(cartj,sphri)
        dimj = nbfxj*ngenj
        call spcart_a_s(buf,scr,dimj,typei,ngeni,FF,printit)
* scr(cartj,sphri) --> buf(sphrj,sphri)
        dimi = nbfsi*ngeni
        call spcart_s_a(scr,buf,dimi,typej,ngenj,FF,printit)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (trani) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartj,carti) --> scr(cartj,sphri)
        dimj = nbfxj*ngenj
        call spcart_a_s(buf,scr,dimj,typei,ngeni,FF,printit)
* scr(cartj,sphri) --> buf(cartj,sphri)
        dimij = nbfxj*ngenj * nbfsi*ngeni
        call dcopy(dimij,scr,1,buf,1)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (tranj) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartj,carti) --> scr(sphrj,carti)
        dimi = nbfxi*ngeni
        call spcart_s_a(buf,scr,dimi,typej,ngenj,FF,printit)
* scr(sphrj,carti) --> buf(sphrj,carti)
        dimij = nbfsj*ngenj * nbfxi*ngeni
        call dcopy(dimij,scr,1,buf,1)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
        write(luout,*)' no basis sets need to be transformed '
        call errquit('spcart_2ctran: fatal error ',911, MEM_ERR)
      endif
*
      end
*.......................................................................
      subroutine spcart_2cBtran(buf,scr,lscr,
     &    nbfxi,nbfsi,typei,ngeni,trani,
     &    nbfxj,nbfsj,typej,ngenj,tranj,
     &    nblocks,printit)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
*
* routine to transform multiple 1e or 2e 2 center blocks of integrals
* buf(jlo:jhi,:ilo:ihi,nblocks) 
* 
      integer lscr           ! [input] length of scratch array
      double precision
     &        scr(lscr)      ! [scratch] scratch array
      integer nbfxi          ! [input] no. of cartesian basis funcs:ish
      integer nbfsi          ! [input] no. of spherical basis funcs:ish
      integer typei          ! [input] angular momentum type:ish
      integer ngeni          ! [input] number general contractions:ish
      integer nbfxj          ! [input] no. of cartesian basis funcs:jsh
      integer nbfsj          ! [input] no. of spherical basis funcs:jsh
      integer typej          ! [input] angular momentum type:jsh
      integer ngenj          ! [input] number general contractions:jsh
      logical trani          ! [input] true -> transform ish block
      logical tranj          ! [input] true -> transform jsh block
      logical printit        ! [input] true -> print transform steps
      integer nblocks        ! [input] number of blocks [intdim,nblock]
      double precision       ! [input/output] cartesian/shperical ints.
     &    buf((nbfxi*ngeni*nbfxj*ngenj),nblocks)
*::local
      integer dimij
      integer dimXX,dimSS
      integer iblock
      logical problem_sp
      logical FF, FT
      FF = .false.
      FT = .true.
      problem_sp =
     &    (typei.eq.-1).and.(ngeni.gt.1)
      problem_sp = problem_sp.and.
     &    (typej.eq.-1).and.(ngenj.gt.1)
      if (problem_sp) then
        write(luout,*)' spcart_2cBtran: sp function error '
        write(luout,*)' sp function block passed with more ',
     &      'than one genereal contraction: ',
     &      'this is not allowed in NWChem'
        write(luout,*)' type i ',typei
        write(luout,*)' ngen i ',ngeni
        write(luout,*)' type j ',typej
        write(luout,*)' ngen j ',ngenj
        call errquit('spcart_2cBtran: fatal error',911, UNKNOWN_ERR)
      endif
*
* sanity checking for spherical transforms
*
      if (trani.and.typei.ge.2.and.nbfxi.le.nbfsi) then
        write(luout,*)' sanity check error on i shell info'
        write(luout,*)' shell type     : ',typei
        write(luout,*)' cartesian size :',nbfxi
        write(luout,*)' spherical size :',nbfsi
        call errquit('spcart_2cBtran:i: fatal error',911, UNKNOWN_ERR)
      endif
      if (tranj.and.typej.ge.2.and.nbfxj.le.nbfsj) then
        write(luout,*)' sanity check error on j shell info'
        write(luout,*)' shell type     : ',typej
        write(luout,*)' cartesian size :',nbfxj
        write(luout,*)' spherical size :',nbfsj
        call errquit('spcart_2cBtran:j: fatal error',911, UNKNOWN_ERR)
      endif
*
* Check for all s, p, sp shells (e.g., no transform required).
*
      dimXX = nbfxi * nbfxj
      dimSS = nbfsi * nbfsj
      if (dimXX.eq.dimSS) return
*
* check scratch size
*
      dimij = nbfxi * ngeni * nbfxj * ngenj
      dimij = dimij * nblocks
      if (dimij.gt.lscr) then
        write(luout,*)' dimij   :',dimij
        write(luout,*)' lscr    :',lscr
        call errquit
     &      ('spcart_2cBtran: error in scratch size for',911, MEM_ERR)
      endif
*
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* transform each block independently
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do iblock = 1 , nblocks
        call spcart_2ctran(buf(1,iblock),scr,lscr,
     &      nbfxi,nbfsi,typei,ngeni,trani,
     &      nbfxj,nbfsj,typej,ngenj,tranj,
     &      printit)
      enddo
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* move each spherical block block down in buf
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dimXX = nbfxi * nbfxj
      dimXX = dimXX * ngeni * ngenj
      dimSS = nbfsi * nbfsj
      dimSS = dimSS * ngeni * ngenj
      call int_c2s_mv(buf,dimXX,dimSS,nblocks,
     &    scr,lscr,'spcart_2cBtran')
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
*.......................................................................
      subroutine spcart_3ctran(buf,scr,lscr,
     &    nbfxi,nbfsi,typei,ngeni,trani,
     &    nbfxj,nbfsj,typej,ngenj,tranj,
     &    nbfxk,nbfsk,typek,ngenk,trank,
     &    printit)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
*::passed
*
* routine to transform a 1e or 2e 3 center block of integrals
* block(klo:khi,jlo:jhi,:ilo:ihi) 
* 
      double precision buf(*)! [input/output] cartesian/shperical ints.
      integer lscr           ! [input] length of scratch array
      double precision
     &        scr(lscr)      ! [scratch] scratch array
      integer nbfxi          ! [input] no. of cartesian basis funcs:ish
      integer nbfsi          ! [input] no. of spherical basis funcs:ish
      integer typei          ! [input] angular momentum type:ish
      integer ngeni          ! [input] number general contractions:ish
      integer nbfxj          ! [input] no. of cartesian basis funcs:jsh
      integer nbfsj          ! [input] no. of spherical basis funcs:jsh
      integer typej          ! [input] angular momentum type:jsh
      integer ngenj          ! [input] number general contractions:jsh
      integer nbfxk          ! [input] no. of cartesian basis funcs:ksh
      integer nbfsk          ! [input] no. of spherical basis funcs:ksh
      integer typek          ! [input] angular momentum type:ksh
      integer ngenk          ! [input] number general contractions:ksh
      logical trani          ! [input] true -> transform ish block
      logical tranj          ! [input] true -> transform jsh block
      logical trank          ! [input] true -> transform ksh block
      logical printit        ! [input] true -> print transform steps
*::local
      integer dimXXX  ! size of full cartesian block
      integer dimSSS  ! size of full spherical block
      integer dimkj, dimji, dimijk  ! intermediate sizes
      integer dimk, dimi            ! intermediate sizes
      logical problem_sp
      logical FF, FT
      FF = .false.
      FT = .true.
      problem_sp =
     &    (typei.eq.-1).and.(ngeni.gt.1)
      problem_sp = problem_sp.and.
     &    (typej.eq.-1).and.(ngenj.gt.1)
      problem_sp = problem_sp.and.
     &    (typek.eq.-1).and.(ngenk.gt.1)
      if (problem_sp) then
        write(luout,*)' spcart_3ctran: sp function error '
        write(luout,*)' sp function block passed with more ',
     &      'than one genereal contraction: ',
     &      'this is not allowed in NWChem'
        write(luout,*)' type i ',typei
        write(luout,*)' ngen i ',ngeni
        write(luout,*)' type j ',typej
        write(luout,*)' ngen j ',ngenj
        write(luout,*)' type k ',typek
        write(luout,*)' ngen k ',ngenk
        call errquit('spcart_3ctran: fatal error',911, UNKNOWN_ERR)
      endif
*
* sanity checking for spherical transforms
*
      if (trani.and.typei.ge.2.and.nbfxi.le.nbfsi) then
        write(luout,*)' sanity check error on i shell info'
        write(luout,*)' shell type     : ',typei
        write(luout,*)' cartesian size :',nbfxi
        write(luout,*)' spherical size :',nbfsi
        call errquit('spcart_3ctran:i: fatal error',911, UNKNOWN_ERR)
      endif
      if (tranj.and.typej.ge.2.and.nbfxj.le.nbfsj) then
        write(luout,*)' sanity check error on j shell info'
        write(luout,*)' shell type     : ',typej
        write(luout,*)' cartesian size :',nbfxj
        write(luout,*)' spherical size :',nbfsj
        call errquit('spcart_3ctran:j: fatal error',911, UNKNOWN_ERR)
      endif
      if (trank.and.typek.ge.2.and.nbfxk.le.nbfsk) then
        write(luout,*)' sanity check error on k shell info'
        write(luout,*)' shell type     : ',typek
        write(luout,*)' cartesian size :',nbfxk
        write(luout,*)' spherical size :',nbfsk
        call errquit('spcart_3ctran:k: fatal error',911, UNKNOWN_ERR)
      endif
*
* Check for all s, p, sp shells (e.g., no transform required).
*
      dimXXX = nbfxi * nbfxj * nbfxk 
      dimSSS = nbfsi * nbfsj * nbfsk 
      if (dimXXX.eq.dimSSS) return
*
* check scratch size
*
      dimijk = nbfxi * ngeni * nbfxj * ngenj * nbfxk * ngenk 
      if (dimijk.gt.lscr) then
        write(luout,*)' dimijk  :',dimijk
        write(luout,*)' lscr    :',lscr
        call errquit
     &      ('spcart_3ctran: error in scratch size for',911,
     &       UNKNOWN_ERR)
      endif
*
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (trani.and.tranj.and.trank) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartk,cartj,carti) --> scr(cartk,cartj,sphri)
        dimkj = nbfxk*ngenk*nbfxj*ngenj
        call spcart_a_s(buf,scr,dimkj,typei,ngeni,FF,printit)
* scr(cartk,cartj,sphri) --> buf(sphrk,cartj,sphri)
        dimji = nbfxj*ngenj*nbfsi*ngeni
        call spcart_s_a(scr,buf,dimji,typek,ngenk,FF,printit)
* buf(sphrk,cartj,sphri) --> scr(sphrk,sphrj,sphri)
        dimk = nbfsk*ngenk
        dimi = nbfsi*ngeni
        call spcart_a_s_b(buf,scr,dimk,dimi,typej,ngenj,FF,printit)
* scr(sphrk,sphrj,sphri) --> buf(sphrk,sphrj,sphri)
        dimijk = nbfsk*ngenk * nbfsj*ngenj * nbfsi*ngeni
        call dcopy(dimijk,scr,1,buf,1)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (trani.and.trank) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartk,cartj,carti) --> scr(cartk,cartj,sphri)
        dimkj = nbfxk*ngenk*nbfxj*ngenj
        call spcart_a_s(buf,scr,dimkj,typei,ngeni,FF,printit)
* scr(cartk,cartj,sphri) --> buf(sphrk,cartj,sphri)
        dimji = nbfxj*ngenj*nbfsi*ngeni
        call spcart_s_a(scr,buf,dimji,typek,ngenk,FF,printit)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (trani.and.tranj) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartk,cartj,carti) --> scr(cartk,cartj,sphri)
        dimkj = nbfxk*ngenk*nbfxj*ngenj
        call spcart_a_s(buf,scr,dimkj,typei,ngeni,FF,printit)
* scr(cartk,cartj,sphri) --> buf(cartk,sphrj,sphri)
        dimk = nbfxk*ngenk
        dimi = nbfsi*ngeni
        call spcart_a_s_b(scr,buf,dimk,dimi,typej,ngenj,FF,printit)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (tranj.and.trank) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartk,cartj,carti) --> scr(sphrk,cartj,carti)
        dimji = nbfxi*ngeni*nbfxj*ngenj
        call spcart_s_a(buf,scr,dimji,typek,ngenk,FF,printit)
* scr(sphrk,cartj,carti) --> buf(sphrk,sphrj,carti)
        dimk = nbfsk*ngenk
        dimi = nbfxi*ngeni
        call spcart_a_s_b(scr,buf,dimk,dimi,typej,ngenj,FF,printit)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (trani) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartk,cartj,carti) --> scr(cartk,cartj,sphri)
        dimkj = nbfxk*ngenk*nbfxj*ngenj
        call spcart_a_s(buf,scr,dimkj,typei,ngeni,FF,printit)
* scr(cartk,cartj,sphri) --> buf(cartk,cartj,sphri)
        dimijk = nbfxk*ngenk * nbfxj*ngenj * nbfsi*ngeni
        call dcopy(dimijk,scr,1,buf,1)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (tranj) then
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* buf(cartk,cartj,carti) --> scr(cartk,sphrj,carti)
        dimk = nbfxk*ngenk
        dimi = nbfxi*ngeni
        call spcart_a_s_b(buf,scr,dimk,dimi,typej,ngenj,FF,printit)
* scr(cartk,sphrj,carti) --> buf(cartk,sphrj,carti)
        dimijk = nbfxk*ngenk * nbfsj*ngenj * nbfxi*ngeni
        call dcopy(dimijk,scr,1,buf,1)
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif (trank) then
* buf(cartk,cartj,carti) --> scr(sphrk,cartj,carti)
        dimji = nbfxi*ngeni*nbfxj*ngenj
        call spcart_s_a(buf,scr,dimji,typek,ngenk,FF,printit)
* scr(sphrk,cartj,carti) --> buf(sphrk,cartj,carti)
        dimijk = nbfsk*ngenk * nbfxj*ngenj * nbfxi*ngeni
        call dcopy(dimijk,scr,1,buf,1)
      else
        write(luout,*)' no basis sets need to be transformed '
        call errquit('spcart_3ctran: fatal error ',911, UNKNOWN_ERR)
      endif
*
      end
*.......................................................................
      subroutine spcart_3cBtran(buf,scr,lscr,
     &    nbfxi,nbfsi,typei,ngeni,trani,
     &    nbfxj,nbfsj,typej,ngenj,tranj,
     &    nbfxk,nbfsk,typek,ngenk,trank,
     &    nblocks,printit)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
*
* routine to transform multiple 1e or 2e 3 center blocks of integrals
* buf(klo:khi,jlo:jhi,:ilo:ihi,nblocks) 
* 
      integer lscr           ! [input] length of scratch array
      double precision
     &        scr(lscr)      ! [scratch] scratch array
      integer nbfxi          ! [input] no. of cartesian basis funcs:ish
      integer nbfsi          ! [input] no. of spherical basis funcs:ish
      integer typei          ! [input] angular momentum type:ish
      integer ngeni          ! [input] number general contractions:ish
      integer nbfxj          ! [input] no. of cartesian basis funcs:jsh
      integer nbfsj          ! [input] no. of spherical basis funcs:jsh
      integer typej          ! [input] angular momentum type:jsh
      integer ngenj          ! [input] number general contractions:jsh
      integer nbfxk          ! [input] no. of cartesian basis funcs:ksh
      integer nbfsk          ! [input] no. of spherical basis funcs:ksh
      integer typek          ! [input] angular momentum type:ksh
      integer ngenk          ! [input] number general contractions:ksh
      logical trani          ! [input] true -> transform ish block
      logical tranj          ! [input] true -> transform jsh block
      logical trank          ! [input] true -> transform ksh block
      logical printit        ! [input] true -> print transform steps
      integer nblocks        ! [input] number of blocks [intdim,nblock]
      double precision       ! [input/output] cartesian/shperical ints.
     &    buf((nbfxi*ngeni*nbfxj*ngenj*nbfxk*ngenk),nblocks)
*::local
      integer dimijk
      integer dimXXX,dimSSS
      integer iblock
      logical problem_sp
      logical FF, FT
      FF = .false.
      FT = .true.
      problem_sp =
     &    (typei.eq.-1).and.(ngeni.gt.1)
      problem_sp = problem_sp.and.
     &    (typej.eq.-1).and.(ngenj.gt.1)
      problem_sp = problem_sp.and.
     &    (typek.eq.-1).and.(ngenk.gt.1)
      if (problem_sp) then
        write(luout,*)' spcart_3cBtran: sp function error '
        write(luout,*)' sp function block passed with more ',
     &      'than one genereal contraction: ',
     &      'this is not allowed in NWChem'
        write(luout,*)' type i ',typei
        write(luout,*)' ngen i ',ngeni
        write(luout,*)' type j ',typej
        write(luout,*)' ngen j ',ngenj
        write(luout,*)' type k ',typek
        write(luout,*)' ngen k ',ngenk
        call errquit('spcart_3cBtran: fatal error',911, UNKNOWN_ERR)
      endif
*
* sanity checking for spherical transforms
*
      if (trani.and.typei.ge.2.and.nbfxi.le.nbfsi) then
        write(luout,*)' sanity check error on i shell info'
        write(luout,*)' shell type     : ',typei
        write(luout,*)' cartesian size :',nbfxi
        write(luout,*)' spherical size :',nbfsi
        call errquit('spcart_3cBtran:i: fatal error',911, UNKNOWN_ERR)
      endif
      if (tranj.and.typej.ge.2.and.nbfxj.le.nbfsj) then
        write(luout,*)' sanity check error on j shell info'
        write(luout,*)' shell type     : ',typej
        write(luout,*)' cartesian size :',nbfxj
        write(luout,*)' spherical size :',nbfsj
        call errquit('spcart_3cBtran:j: fatal error',911, UNKNOWN_ERR)
      endif
      if (trank.and.typek.ge.2.and.nbfxk.le.nbfsk) then
        write(luout,*)' sanity check error on k shell info'
        write(luout,*)' shell type     : ',typek
        write(luout,*)' cartesian size :',nbfxk
        write(luout,*)' spherical size :',nbfsk
        call errquit('spcart_3cBtran:k: fatal error',911, UNKNOWN_ERR)
      endif
*
* Check for all s, p, sp shells (e.g., no transform required).
*
      dimXXX = nbfxi * nbfxj * nbfxk 
      dimSSS = nbfsi * nbfsj * nbfsk 
      if (dimXXX.eq.dimSSS) return
*
* check scratch size
*
      dimijk = nbfxi * ngeni * nbfxj * ngenj * nbfxk * ngenk 
      dimijk = dimijk * nblocks
      if (dimijk.gt.lscr) then
        write(luout,*)' dimijk  :',dimijk
        write(luout,*)' lscr    :',lscr
        call errquit
     &      ('spcart_3cBtran: error in scratch size for',911,
     &       UNKNOWN_ERR)
      endif
*
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* transform each block independently
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do iblock = 1 , nblocks
        call spcart_3ctran(buf(1,iblock),scr,lscr,
     &      nbfxi,nbfsi,typei,ngeni,trani,
     &      nbfxj,nbfsj,typej,ngenj,tranj,
     &      nbfxk,nbfsk,typek,ngenk,trank,
     &      printit)
      enddo
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* move each spherical block block down in buf
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dimXXX = nbfxi * nbfxj * nbfxk
      dimXXX = dimXXX * ngeni * ngenj * ngenk
      dimSSS = nbfsi * nbfsj * nbfsk
      dimSSS = dimSSS * ngeni * ngenj * ngenk
      call int_c2s_mv(buf,dimXXX,dimSSS,nblocks,
     &    scr,lscr,'spcart_3cBtran')
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
