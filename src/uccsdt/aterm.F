      subroutine uccsdt_aterm(urange,vrange)
c
c$Id: aterm.F,v 1.4 2002-10-05 01:22:26 bert Exp $
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
      integer list_offset,mnpairs,m,n
      integer n_address,m_address
      integer g_t2_iejf, g_x
      integer l_vulist,k_vulist,l_upart,k_upart,l_vpart,k_vpart
      double precision prefactor
C
C     GENERATION OF PART OF A-TERM, STORED IN X(1-3) ARRAYS
C
c  R(me,nf) aa = 0.5* t(ie,jf)[<ij|mn>+0.5*X1(ij,mn) i=j=m=n=e=f=alpha
c  R(me,nf) bb = 0.5* t(ie,jf)[<ij|mn>+0.5*X2(ij,mn) i=j=m=n=e=f=beta
c  R(me,nf) ab =      t(ie,jf)[<ij|mn>+    X3(ij,mn) m=e=i=alpha n=f=j=beta
c
c  X1(ij,mn) = <ij||ab>*t(ma,nb) -> can be generated from class 1 ??? i=j=a=b=m=n=alpha
c  X2(ij,mn) = <ij||ab>*t(ma,nb) -> can be generated from class 4 ??? i=j=a=b=m=n=beta
c  X3(ij,mn) = <ij||ab>*t(ma,nb) -> can be generated from class 2 ??? i=m=a=alpha j=n=b=beta
c
c
c  We can do them by order of R or by order of X, by R is less memory (only one X memory)
c
      do nspin = 1, 2
         do mspin = 1, nspin
            list_offset = nspin + (mspin-1)*2
            if (nspin .eq. mspin) then
               prefactor = 0.5d0
            else
               prefactor = 1.0d0
            do symnj = 0, nir-1 
               symnv = symnj
               symmi = symnj 
               symmu = symmi
               ispin = mspin
               jspin = nspin
               dimmi = ov_len(symmi,mspin,ispin)
               dimnj = ov_len(symnj,nspin,jspin)
               if (.not. ga_create(MT_DBL,dimmi,dimnj,'X_temp',0,0,
     &             g_x) call errquit()
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
               call ga_scale(g_x,prefactor)
c
c     Add <ij|mn> integrals from classes 5,6,8
c
c     here use list_offset+4
               list_offset = list_offset + 4
               do symj = 0, nir-1
                  symv = symj
                  symn = ieor(symnj, symj)
                  do symi = 0, nir-1
                     symu = symi
                     symm = ieor (symmi, symi)
                     dimi = nv_sym(symi,ispin)
                     dimj = nv_sym(symj,jspin)
                     dimm = nv_sym(symm,mspin)
                     dimn = nv_sym(symn,nspin)
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
     &                      'vu list',l_vulist,k_vulist))
                        dim_mx = max(dimi,dimj)
                        if (.not. ma_push_get(mt_dbl,dim_mx*dim_mx,
     &                              'ij list',l_ijlist, k_ijlist))
                        n_address = 1+n-v_sym(1,symn,spinn)
     &                              +oso_off(v,symnv,spinn)
                        m_address = 1+m-v_sym(1,symm,spinm)
     &                              +oso_off(u,symmu,spinm)
                        call ga_get(list(symmu,list_nr),
     &                              n_address,dimv,nv_sym(symn,1),
     &                              m_address,dimu,nv_sym(symm,1),
     &                              dbl_mb(k_vulist),dimv)
                        if (.not. ma_push_get(mt_dbl,dimu*dimi,
     &                      'c part ui',l_upart,k_upart))
                        if (.not. ma_push_get(mt_dbl,dimv*dimj,
     &                      'c part vj',l_vpart,k_vpart))
                        call ga_get(g_part(ispin),ulo,uhi,
     &                      v_sym(1,symi,ispin),v_sym(2,symi,espin),
     &                      k_upart,ld)
                        call ga_get(g_part(jspin),vlo,vhi,
     &                      v_sym(1,symj,jspin),v_sym(2,symj,fspin),
     &                      k_vpart,ld)
                        call dgemm('n','n',dimv,dimi,dimu,1.0d0,
     &                      dbl_mb(k_vulist),dimv,dbl_mb(k_upart),dimu,
     &                      1.0d0,dbl_mb(k_ijlist),dimv)
                        call dgemm('t','n',dime,dimj,dimv,1.0d0,
     &                      dbl_mb(k_ijlist),dimv,dbl_mb(K_vpart),dimv,
     &                      0.0d0,dbl_mb(k_ijlist),dimi)
                        if (.not. ma_pop_stack(l_upart)) call 
     &                      errquit()
                        if (.not. ma_pop_stack(l_vpart)) call 
     &                      errquit()
                        call ga_put(g_x,n_address,dimi,dimm,
     &                      m_address,dimj,dimn,dbl(k_ijlist),dimi,
     &                      1.0d0)
                        if (.not. ma_pop_stack(l_vulist)) call 
     &                      errquit()
                        if (.not. ma_pop_stack(l_ijlist)) call 
     &                      errquit()
                     end do
                  end do
               end do
               call get_T(g_t2_iejf,symnj,ispin,espin,jspin,fspin)
c              prefactor*t(ie,jf)*x(ij,mn) -> r(me,nf)
            end do
         end do
      enddo
c
      end

