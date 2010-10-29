      subroutine uccsdt_bgxterms(d_t,d_b,d_g,d_x)
c $Id$
c
c     Compute B, X, and G terms defined as follows:
c
c     Pure spin
c
c     B(a,b,i,j) = 1/2 <ab|^|cd> t(c,d,i,j)
c
c     ... t(u,v,i,j) = C(u,c)C(v,d)t(c,d,i,j)
c
c     ... Omega(u,v,i,j) = {<uv||wl> t(w,l,i,j)}/(1+delta_wl) ; w>=l
c
c     ... B(a,b,i,j) = Omega(u,v,i,j)C(u,a)C(v,b)
c
c     X(i,j,m,n) = <ij|^|ab> t(a,b,m,n)
c
c     ... X(i,j,m,n) = Omega(u,v,m,n)C(u,i)C(v,j)
c
c     G(e,m) = 1/2 <ie|^|ab> t(e,b,i,j)
c                    _    _      _   _
c                + <ei|^ ab> t(a,b,m,i)
c
c     ... G(e,m) = Omega(u,v,i,m)C(u,e)C(v,i) 
c                              _           _
c     ...        + Omega(u,v,m,i)C(u,e)C(v,i)
c
c     Mixed spin
c         _   _      _    _      _   _
c     B(a,b,i,j) = <ab|^ cd> t(c,d,i,j)
c                 _              _     _   _
c     ... t(u,v,i,j) = C(u,c)C(v,d)t(c,d,i,j)
c                      _            _            _
c     ... t[+/-](u,v,i,j) = t(u,v,i,j)+/-t(v,u,i,j)
c                     _                            _
c     ... Omega(u,v,i,j) = 1/2 {<uv||wl>[+]t(w,l,i,j)[+]
c                                                  _
c     ...                     + <uv||wl>[-]t(w,l,i,j)[-]}/(1+delta_wl) ; w>=l
c             _   _                _           _
c     ... B(a,b,i,j) = Omega(u,v,i,j)C(u,a)C(v,b)
c         _   _      _    _      _   _
c     X(i,j,m,n) = <ij|^ ab> t(a,b,m,n)
c             _   _                _           _
c     ... X(i,j,m,n) = Omega(u,v,m,n)C(u,i)C(v,j)
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "global.fh"
      integer d_t           ! [input]  DRA handle for t(a,b,i,j)
      integer d_b           ! [output] DRA handle for b(a,b,i,j)
      integer d_x           ! [output] DRA handle for x(a,b,i,j)
      integer d_g           ! [output] DRA handle for g(a,i)
      integer g_uv_p        ! GA handle for t[+](u,v,i,j)
      integer g_uv_m        ! GA handle for t[-](u,v,i,j)
      integer g_omega       ! GA handle for Omega(u,v,i,j)
      integer spini,spinj
c
c     Spin loop
c
csh   if (.not. uccsdt_ampfile_zero(d_g)) 
csh  $ call errquit('uccsdt_bxgterms: ampfile_zero?',999)
      do spini = 0,1
       do spinj = spini,1
c
c     SO half transform
c
        call uccsdt_tosopair(d_t,g_uv_p,g_uv_m,spini,spinj)
c
c     Contract
c
        call uccsdt_b_contract(g_uv_p,g_uv_m,g_omega,spini,spinj)
c
c     Make B
c
        call uccsdt_makeb(g_omega,d_b,spini,spinj)
c
c     Make X
c
        call uccsdt_makex(g_omega,d_x,spini,spinj)
c
c     Make G
c
        call uccsdt_makeg(g_omega,d_g,spini,spinj)
c
       end do
      end do
      return
      end
