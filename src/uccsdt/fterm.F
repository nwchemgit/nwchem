      subroutine uccsdt_fterm(urange,vrange,lists)
c
c$Id$
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7),vrange(2,0:7),u,v,m,n,vhi,vlo,uhi,ulo
      integer nspin,mspin,espin,fspin
      integer symv,symn,symf,symm,syme,symu
      integer symiu,symmu,symnf,symnv,symme
      integer k_vulist,l_vulist
      integer k_upart,l_upart
      integer k_vpart,l_vpart
      integer k_fpart,l_fpart
      integer k_eflist,l_eflist
      integer m_address,n_address, offset_me,offset_n,offset_nf
      integer dimu,dimv,dimf,dimn,dimm,dime,dimme,dimnf
      integer g_r_array,g_r
      integer dim_mx,list_nr,lists(0:7,14),mnpairs,mn
C
C     GENERATION OF F-TERM
C
C     COMPUTING F(me,nf) = <ef||^ mn> for each spin-block
C
C     <uv|mn> stored as (nv,mu) to be transformed to (nf,me)
C     On each processor
C        Get an mn pair
C        For each mn pair, get whole (u,v)_mn block and transform to (e,f)_mn
c            (v,u)*c(u,e)->(v,e)->(v,e)^t*(v,f)->(e,f)->(me,nf)
C        Put (e,f)_mn block back into F(me,nf)
C
      do nspin = 1, 2
         do mspin = 1, nspin
            espin = mspin
            fspin = nspin
            list_nr = 4 + nspin + (mspin-1)*2
            do symnf = 0, nir-1 
               symnv = symnf
               symme = symnf 
               symmu = symme
               dimme = ov_len(symme,mspin,mspin)
               dimnf = ov_len(symnf,nspin,nspin)
               if (.not. ga_create(MT_DBL,dimme,dimnf,'R_temp',0,0,
     &             g_r_array)) call 
     &             errquit('fterm: create g_r_array failed',0)
               do symf = 0, nir-1
                  symv = symf
                  symn = ieor(symnf, symf)
                  do syme = 0, nir-1
                     symu = syme
                     symm = ieor (symme, syme)
                     dime = nv_sym(syme,espin)
                     dimf = nv_sym(symf,fspin)
                     dimm = no_sym(symm,mspin)
                     dimn = no_sym(symn,nspin)
                     ulo = urange(1,symu)
                     uhi = urange(2,symu)
                     vlo = vrange(1,symv)
                     vhi = vrange(2,symv)
                     dimv = vhi - vlo + 1
                     dimu = uhi - ulo + 1
                     mnpairs = dimm*dimn
                     do mn = ga_nodeid()+1, mnpairs, ga_nnodes() 
                        m = (mnpairs / (dimn+1)) + 1
                        n = mnpairs - (dimm * (m-1))
                        if (.not. ma_push_get(mt_dbl,dimv*dimu,
     &                      'vu list',l_vulist,k_vulist)) call
     &                      errquit('fterm: vu list alloc failed',0)
                        dim_mx = max(dime,dimf)
                        if (.not. ma_push_get(mt_dbl,dim_mx*dim_mx,
     &                      'ef list',l_eflist, k_eflist)) call
     &                      errquit('fterm: ef list alloc failed',0)
                        n_address = o_sym(1,symn,nspin)+n-1
     &                              +oso_v_off(v,symnv,nspin)
                        m_address = o_sym(1,symm,mspin)+m-1
     &                              +oso_u_off(u,symmu,mspin)
                        call ga_get_strided(lists(symmu,list_nr),
     &                             n_address,dimv,dimn,
     &                             m_address,dimu,dimm,
     &                             dbl_mb(k_vulist))
                        if (.not. ma_push_get(mt_dbl,dimu*dime,
     &                      'c particle ue',l_upart,k_upart)) call
     &                      errquit('fterm: ue part alloc failed',0)
                        if (.not. ma_push_get(mt_dbl,dimv*dimf,
     &                      'c particle vf',l_vpart,k_vpart)) call
     &                      errquit('fterm: vf part alloc failed',0)
                        call ga_get(g_part(espin),ulo,uhi,
     &                             v_sym(1,syme,espin),
     &                             v_sym(2,syme,espin),
     &                             dbl_mb(k_upart),dimu)  
                        call ga_get(g_part(fspin),vlo,vhi,
     &                             v_sym(1,symf,fspin),
     &                             v_sym(2,symf,fspin),
     &                             dbl_mb(k_vpart),dimv)  
                        call dgemm('n','n',dimv,dime,dimu,1.0d0,
     &                             dbl_mb(k_vulist),dimv,
     &                             dbl_mb(k_upart),dimu,1.0d0,
     &                             dbl_mb(k_eflist),dimv)
                        call dgemm('t','n',dime,dimf,dimv,1.0d0,
     &                             dbl_mb(k_eflist),dimv,
     &                             dbl_mb(k_vpart),dimv,1.0d0,
     &                             dbl_mb(k_eflist),dime)
                        if (.not. ma_pop_stack(l_vpart)) call 
     &                      errquit('fterm: vf part dealloc failed',0)
                        if (.not. ma_pop_stack(l_upart)) call 
     &                      errquit('fterm: ue part dealloc failed',0)
                        call ga_put_strided(g_r_array,
     &                             m_address,dime,dimm,
     &                             n_address,dimf,dimn,
     &                             dbl_mb(k_eflist),dime)
                        if (.not. ma_pop_stack(l_eflist)) call 
     &                      errquit('fterm: ef list dealloc failed',0)
                        if (.not. ma_pop_stack(l_vulist)) call 
     &                      errquit('fterm: vu list dealloc failed',0)
                     end do
                  end do
               enddo
               offset_me = file_t2_off(mspin,espin,nspin,fspin,symnf)+
     &                     ov_off(1,symme,mspin,espin)
               offset_nf = file_t2_off(mspin,espin,nspin,fspin,symnf)+
     &                     ov_off(1,symnf,nspin,fspin)
               call ga_copy_patch('n',g_r_array,1,dimme,1,dimnf,g_r,
     &              offset_me,offset_me+dimme,offset_nf,offset_nf+dimnf)
               if (.not. ga_destroy(g_r_array)) call
     &             errquit('fterm: g_r_array dealloc failed',0)
            end do
         end do
      enddo
c
      end
c
      subroutine ga_get_strided(ga_in,ilo,ilength,istride,
     &                          jlo,jlength,jstride,array_out)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
      integer ga_in,ilo,ilength,istride,jlo,jlength,jstride
      integer nelem,ilist,jlist,ihi,jhi
      integer k_i,l_i,k_j,l_j
      double precision array_out(*)
c
      ihi = ilo + (ilength-1) * istride
      jhi = jlo + (jlength-1) * jstride
      nelem = ilength*jlength
      if (.not. ma_push_get(mt_dbl,nelem,'i list',l_i,k_i)) call
     &    errquit('ga_get_strided: i list alloc failed',0)
      if (.not. ma_push_get(mt_dbl,nelem,'j list',l_j,k_j)) call
     &    errquit('ga_get_strided: j list alloc failed',0)
      nelem = 0
      do jlist = jlo, jhi, jstride
         do ilist = ilo, ihi, istride
            int_mb(k_i+nelem) = ilist 
            int_mb(k_j+nelem) = jlist 
            nelem = nelem + 1
         enddo
      enddo
      call ga_gather(ga_in,array_out,int_mb(k_i),int_mb(k_j),nelem)
      if (.not. ma_pop_stack(l_i)) call 
     &   errquit('ga_put_strided: i list dealloc failed',0)
      if (.not. ma_pop_stack(l_j)) call 
     &   errquit('ga_put_strided: j list dealloc failed',0)
c
      end
c
      subroutine ga_put_strided(ga_out,ilo,ilength,istride,
     &                          jlo,jlength,jstride,array_in)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
      integer ga_out,ilo,ilength,istride,jlo,jlength,jstride
      integer nelem,ilist,jlist,ihi,jhi
      integer k_i,l_i,k_j,l_j
      double precision array_in(*)
c
      ihi = ilo + (ilength-1) * istride
      jhi = jlo + (jlength-1) * jstride
      nelem = ilength*jlength
      if (.not. ma_push_get(mt_dbl,nelem,'i list',l_i,k_i)) call
     &    errquit('ga_put_strided: i list alloc failed',0)
      if (.not. ma_push_get(mt_dbl,nelem,'j list',l_j,k_j)) call
     &    errquit('ga_put_strided: j list alloc failed',0)
      nelem = 0
      do jlist = jlo, jhi, jstride
         do ilist = ilo, ihi, istride
            nelem = nelem + 1
            int_mb(k_i+nelem) = ilist 
            int_mb(k_j+nelem) = jlist 
         enddo
      enddo
      call ga_scatter(ga_out,array_in,int_mb(k_i),int_mb(k_j),nelem)
      if (.not. ma_pop_stack(l_i)) call 
     &   errquit('ga_put_strided: i list dealloc failed',0)
      if (.not. ma_pop_stack(l_j)) call 
     &   errquit('ga_put_strided: j list dealloc failed',0)
c
      end
