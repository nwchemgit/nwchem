      subroutine uccsdt_aterm(urange,vrange)
c
c$Id: aterm.F,v 1.7 2002-10-22 23:51:17 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7), vrange(2,0:7)
      integer nspin,mspin,ispin,jspin
      integer symnj,symnv,symmi,symmu,symj,symv,symn,symi,symu,symm
      integer dimmi,dimnj,dimi,dimj,dimu,dimv,dimm,dimn
      integer ulo,vlo
      integer list_offset,list_nr,mnpairs,m,n
      integer n_address,m_address
      integer g_t2_iejf, g_x
      integer l_vulist,k_vulist,l_upart,k_upart,l_vpart,k_vpart
      double precision prefactor
C
C     GENERATION OF PART OF A-TERM, STORED IN X(1-3) ARRAYS
C
c  R(me,nf) aa = 0.5* t(ie,jf)[<ij|mn>+0.5*X1(ij,mn)] i=j=m=n=e=f=alpha
c  R(me,nf) bb = 0.5* t(ie,jf)[<ij|mn>+0.5*X2(ij,mn)] i=j=m=n=e=f=beta
c  R(me,nf) ab =      t(ie,jf)[<ij|mn>+    X3(ij,mn)] m=e=i=alpha n=f=j=beta
c
c  X1(ij,mn) = <ij||ab>*t(ma,nb) -> can be generated from class 1 ??? i=j=a=b=m=n=alpha
c  X2(ij,mn) = <ij||ab>*t(ma,nb) -> can be generated from class 4 ??? i=j=a=b=m=n=beta
c  X3(ij,mn) = <ij||ab>*t(ma,nb) -> can be generated from class 2 ??? i=m=a=alpha j=n=b=beta
c
c  So will have X-terms in 3 GA's of format (m,i,symi,n,j,symj,nj) for all nj's and I'll use those. 
c  Need to get the offsets from So's Xterm code.
c     Scale X1 and X2
c     Add <ij|mn> integrals in correct places
c     Get complete t(ie,jf) spin-block. How to matmul those? Probably like: 
c        collect (i,j)_mn from X
c        collect (i,j)_ef from T2
c        matmul blocks into R(
c     
c  We can do them by order of R or by order of X, by R is less memory (only one X memory)
c
      do nspin = 1, 2
         do mspin = 1, nspin
            list_offset = nspin + (mspin-1)*2
            if (nspin .eq. mspin) then
               prefactor = 0.5d0
               g_r = g_rbb
               if (nspin .eq. 1) g_r = g_raa
            else
               prefactor = 1.0d0
               g_r = g_rab
            do symnf = 0, nir-1 
               symme = symnf
               symmi = symnf 
               symmu = symmi
               ispin = mspin
               espin = mspin
               jspin = nspin
               fspin = nspin
               dimme = ov_len(symme,mspin,espin)
               dimnf = ov_len(symnf,nspin,fspin)
c
c     Generate X1, X2 and X3 from classes 1,2,4 -> probbaly cannot be done as those 
c     classes are <uv||ij> and we need <ij||uv> !!!!!
c
c     here use list_offset
c
c     So has already made the X-terms in his B-term routines. Maybe I should do a partial
c     A term here while so does the other part coming from the X-term, read them in from 
c     GA or disk
c
               call ga_scale(g_x(list_offset),prefactor)
c
c     Add <ij|mn> integrals from classes 5,6,8
c
c     here use list_offset+4
               list_nr = list_offset + 4
               do symj = 0, nir-1
                  symv = symj
                  symn = ieor(symnj, symj)
                  do symi = 0, nir-1
                     symu = symi
                     symm = ieor (symmi, symi)
                     dimi = no_sym(symi,ispin)
                     dimj = no_sym(symj,jspin)
                     dimm = no_sym(symm,mspin)
                     dimn = no_sym(symn,nspin)
                     ulo = urange(1,symu)
                     uhi = urange(2,symu)
                     vlo = vrange(1,symv)
                     vhi = vrange(2,symv)
                     dimv = vhi - vlo + 1
                     dimu = uhi - ulo + 1
                     mnpairs = no_sym(symm,mspin)*no_sym(symn,nspin)
                     do mn = ga_nodeid()+1, mnpairs, ga_nnodes()
                        m = (mnpairs / (no_sym(symn,nspin)+1)) + 1
                        n = mnpairs - (no_sym(symm,mspin) * (m-1))
                        if (.not. ma_push_get(mt_dbl,dimv*dimu,
     &                      'vu list',l_vulist,k_vulist)) call
     &                      errquit('aterm: k_vulist alloc failed',0)
                        dim_mx = max(dimi,dimj)
                        if (.not. ma_push_get(mt_dbl,dim_mx*dim_mx,
     &                      'ij list',l_ijlist, k_ijlist)) call
     &                      errquit('aterm: k_ijlist alloc failed',0)
                        n_address = o_sym(1,symn,nspin)+n-1
     &                              +oso_v_off(v,symnv,nspin)
                        m_address = o_sym(1,symm,mspin)+m-1
     &                              +oso_u_off(u,symmu,mspin)
                        call ga_get(list(symmu,list_nr),
     &                              n_address,dimv,no_sym(symn,1),
     &                              m_address,dimu,no_sym(symm,1),
     &                              dbl_mb(k_vulist),dimv)
                        if (.not. ma_push_get(mt_dbl,dimu*dimi,
     &                      'c part ui',l_upart,k_upart)) call
     &                      errquit('aterm: k_upart alloc failed',0)
                        if (.not. ma_push_get(mt_dbl,dimv*dimj,
     &                      'c part vj',l_vpart,k_vpart)) call
     &                      errquit('aterm: k_vpart alloc failed',0)
                        call ga_get(g_part(ispin),ulo,uhi,
     &                      o_sym(1,symi,ispin),o_sym(2,symi,ispin),
     &                      k_upart,ld)
                        call ga_get(g_part(jspin),vlo,vhi,
     &                      o_sym(1,symj,jspin),o_sym(2,symj,jspin),
     &                      k_vpart,ld)
                        call dgemm('n','n',dimv,dimi,dimu,1.0d0,
     &                      dbl_mb(k_vulist),dimv,dbl_mb(k_upart),dimu,
     &                      1.0d0,dbl_mb(k_ijlist),dimv)
                        call dgemm('t','n',dime,dimj,dimv,1.0d0,
     &                      dbl_mb(k_ijlist),dimv,dbl_mb(k_vpart),dimv,
     &                      0.0d0,dbl_mb(k_ijlist),dimi)
                        if (.not. ma_pop_stack(l_vpart)) call 
     &                      errquit('aterm: k_vpart dealloc failed',0)
                        if (.not. ma_pop_stack(l_upart)) call 
     &                      errquit('aterm: k_upart dealloc failed',0)
C
C    CANNOT DO A PUT AS G_X ALREADY HAS X-TERM IN IT. 
C    1. DO A GA_ACC (BUT IS THIS STRIDED)
C    2. DO A GA_GET_STRIDED FROM G_X, ADD NEW CONTRIBUTION AND DO GA_PUT_STRIDED
C
                        call ga_put(g_x(list_offset),n_address,dimi,
     &                      dimm,m_address,dimj,dimn,dbl(k_ijlist),
     &                      dimi,1.0d0)
                        if (.not. ma_pop_stack(l_vulist)) call 
     &                      errquit('aterm: k_vulist dealloc failed',0)
                        if (.not. ma_pop_stack(l_ijlist)) call 
     &                      errquit('aterm: k_ijlist dealloc failed',0)
C23456789012345678901234567890123456789012345678901234567890123456789012
                     end do
                  end do
               end do
c              prefactor*t(ie,jf)*X(nj,mi) -> r(me,nf)
c              need only occupied f's
               flo = 
               fhi =
               call uccsdt_ampfile_read_t2(D_AMP,ispin,espin,jspin,
     &                                     fspin,symnj,flo,fhi,
     &                                     g_t2_iejf,.true.,'block')
               call ga_dgemm('n','n',dimnj,dimnj,dimnj,prefactor,
     &                       g_t2_iejf,g_x,1.0d0,g_r)
               if (.not. ga_destroy(g_t2_iejf)) call
     &             errquit('aterm: g_t2_iejf dealloc failed',0)
               if (.not. ga_destroy(g_x))) call
     &             errquit('aterm: g_x dealloc failed',0)
            end do
         end do
      enddo
c
      end

