      subroutine uccsdt_fterm(urange,vrange)
c
c$Id: fterm.F,v 1.1 2002-08-28 22:18:19 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7), vrange(2,0:7)
c
c     i,j = particle-transformed occupied orbitals
c     m,n = hole-transformed occupied orbitals
c     a,b,e,f = holes
c     u,v = symmetry-adapted occupied orbitals
c
c     u - the SOs arising from the unique AO shells in the
c     .   range ushuqlo, ushuqhi
c     All other indices span the complete range.
c
c     IN CCTRANS:
c     #  Integral        Spins           Storage  Used for
c     -- --------------- --------------- -------- --------------
c     1. (ui|vj)-(uj|vi) i=alpha j=alpha I(jv,iu) Z1,Z6
c     2. (ui|vj)         i=alpha j=beta           Z2,Z7
c     3. (ui|vj)         i=beta  j=alpha          Z5,Z8
c     4. (ui|vj)-(uj|vi) i=beta  j=beta           Z3,Z4
c
c     IN CTERM:
c     Lists 1,2,3,4 are stored as the amplitudes but as if u (SO)
c     were a virtual orbital in a reduced space.  And each symmetry
c     block of jb needs to be in a distinct GA (or copied there when
c     used).
c     
c     (note v->b must be done before fn loop)
c
c     1. (iu|jb) i=a=alpha j=b=alpha   antisymmetrized
c     2. (iu|jb) i=a=alpha j=b=beta    
c     3. (iu|jb) i=a=beta  j=b=alpha
c     4. (iu|jb) i=a=beta  j=b=beta    antisymmetrized
c
c
c     IN CCTRANS:
c     5. (um|vn)-(un|vm) m=alpha n=alpha I(nv,mu) aaA and F
c     6. (um|vn)         m=alpha n=beta           abA and F
c    *7. (um|vn)         m=beta  n=alpha          *Not used
c     8. (um|vn)-(un|vm) m=beta  n=beta           bbA and F
c
c     IN CTERM:
c     Lists 5,6,7,8 stored as <uv|mn> distributed over mn.
c
c     5. (mu|nf) m=alpha n=alpha       antisymmetrized
c     6. (mu|nf) m=alpha n=beta
c     7. (mu|nf) m=beta  n=alpha ... NOT NEEDED?
c     8. (mu|nf) m=beta  n=beta        antisymmetrized
c
c
c     IN CCTRANS:
c     9. (ui|vn)-(uv|in) i=alpha n=alpha I(nv,iu) aaC1,abC4
c     10.(ui|vn)         i=alpha n=beta           abC1,bbC2
c     11.(ui|vn)         i=beta  n=alpha          abC5,aaC2
c     12.(ui|vn)-(uv|in) i=beta  n=beta           abC1,abC2
c
c     13.(uv|in)         i=alpha n=alpha I(nv,iu) abC3
c     14.(uv|in)         i=beta  n=beta           abC5
c
c     IN CTERM:
c     Lists 9,10,11,12,13,14 also stored as the amplitudes.
c     The antisymmetrization must take place before the
c     final index transformation of lists 13 and 14.
c
c     9. I(iu,nf) = <if||un> = (iu|nf) - (in|uf) i=n=alpha n=f=alpha 
c     10.I(iu|nf) = <if|un>  = (iu|nf) i=u=alpha n=f=beta 
c     11.I(iu|nf) = <if|un>  = (iu|nf) i=u=beta  n=f=alpha 
c     12.I(iu,nf) = <if||un> = (iu|nf) - (in|uf) i=beta  n=beta      
c     13.(in|uv) i=n=alpha 
c     14.I(iu,nf) = <fi|un> = (in|uf) i=n=beta f=u=alpha
c
c     IN CCTRANS:
c     The occ-SO pairs are stored as follows
c
c     fill oso_off with -99999999
c     ind = 0
c     do symiu
c     .  do u in natural order
c     .     -> symu and symi
c     .     oso_off(u,symiu) = ind
c     .     do i of symi
c     .        pair(1+ind) = T(i,u)
c     .        ind = ind + 1
c    
c     Can also address pairs as
c     pair(1 + i-o_sym(1,symi,spini) + oso_off(u,symiu)) = T(i,u)
c
c     The lists will be used in large matrix multiplications
c     as follows:
c
c     Z(iu,nw) = <ij|uv>*T(jv,nw)
c     do sym(iu) -> sym(jv) -> sym(nw)
c     .   Read T, allocate Z
c     .   Z(iu,nw) <- I(jv,iu)*T(jv,nw)
c     end do
c
c     For best performance we need the symmetry sub-blocks distributed
c     across the whole machine as separate dense arrays.  Thus, each
c     list is stored with a separate GA for each pair symmetry.
c

c
c     STORED LIKE T(nv,mu) -> (um||vn) -> <uv||mn>
c     Have lists 5-8 <uv|mn> distributed over mn ... transform
c     uv to ef, and accumlate <ef|mn> to disk.
c     <ef||mn> -> R(me,nf)
c

C
C     GENERATION OF F-TERM
C
      do nspin = 1, 2
         do mspin = 1, nspin
c           cases do aa (list 5), bb (list 8) and ab (list 6)
            list_nr = 4 + nspin + (mspin-1)*2
            do symnf = 0, nir-1 
               symnv = symnf
               symme = symnf 
               symmu = symme
c              create a temp result vector for R(me,nf) (aa,aa)
c              menf = aaaa, aabb or bbbb, hence mspin=espin, nspin=fspin
               dimme = ov_len(symme,mspin,mspin)
               dimnf = ov_len(symnf,nspin,nspin)
               if (.not. ga_create(MT_DBL,dimme,dimnf,'R_temp',0,0,
     &             g_r_array) call errquit()
               do symf = 0, nir-1
                  symv = symf
                  symn = ieor(symnf, symf)
c                 -> sym(n) -> sym(m)
c     want to transform <uv|mn> to <ef|mn>=R(me,nf) -> R<fe|nm> = R<ef|mn>
c
c     Assume for now that can hold all amplitudes for
c     given symmetry subblock????????
c
c     symv = symf and symu = syme
c
c                 do blocks of f (hold R(me,nf) in memory)
c                  do f = flo(symf,iset), fhi(symf,iset), f_block_size
c                  NO BLOCKS, we have all of T and R in memory
c                  R(m,e,syme,n,f,symf,symnf=symme)
c                  better R(n,f,symf,m,e,syme,symme=symnf)
c
c                    allocate R 
c                    zero or get R, get because ublock gives only part of e ????
c
c                    do mu in parallel mu for given sym(u) and sym(m)
c                       symmu = symme = symnf = symnv
c                       do symu = 0, nir-1
c                          -> symm = ieor(symme, symu)
c
c                    integrals stored as I(nv,mu) -> I(nf,me)
c nv, mu -> n and m holes
c  
                  do syme = 0, nir-1
                     symu = syme
                     symm = ieor (symme, syme)
                     mnpairs = no_sym(symm,mspin)*no_sym(symn,nspin)
                     do mn = 1, mnpairs, mynodeid+1
                        dimv = vrange(2,symv)-vrange(1,symv)+1
                        dimu = urange(2,symu)-urange(1,symu)+1
c                       spin cases e (a in aaaa case and abab case, b in bbbb case)
c                                  f (a in aaaa case, b in abab case and bbbb case)
c                                  e has same spin as m, f has same spin as n
                        espin = mspin
                        fspin = nspin
                        dime = nv_sym(syme,espin)
                        dimf = nv_sym(symf,fspin)
                        dimm = nv_sym(symm,mspin)
                        dimn = nv_sym(symn,nspin)
                        lennv_a = iso_len(symnv,0)
                        lennv_b = iso_len(symnv,1)
                        leniu_a = iu_len(symiu,0)
                        leniu_b = iu_len(symiu,1)
                        ulo = urange(1,syme)
                        uhi = urange(2,syme)
                        vlo = vrange(1,symv)
                        vhi = vrange(2,symv)
c                       allocate ma with correct dimensions
                        if (.not. ma_push_get(mt_dbl,dimv*dimu,
     &                             'vu list',l_vulist,k_vulist))
                        dim_mx = max(dime,dimf)
                        if (.not. ma_push_get(mt_dbl,dim_mx*dim_mx,
     &                             'ef list',l_eflist, k_eflist))
c                       get address (n,1,m,u) for all v -> (v,u)
                        n_address = 1+n-v_sym(1,symn,spinn)
     &                              +oso_off(v,symnv,spinn)
                        m_address = 1+m-v_sym(1,symm,spinm)
     &                              +oso_off(u,symmu,spinm)
                        call ga_get(list(symmu,list_nr),
     &                             n_address,dimv,dimn,
     &                             m_address,dimu,dimm,
     &                             dbl_mb(k_vulist),dimv)
c                    (n,v,m,u) -> (n,f,m,e)
c                    (v,u)*c(u,e)->(v,e)->(v,e)^t*(v,f)->(e,f)->(me,nf)
c                    now you have all uv local for one nm
c                    transform uv to ef and put it into R
c                    
c                    get data from g_hole into MA and do a dgemm
c                    first urange,all_e
c                    second vrange, all_f
c                    same for next operation for <ij||mn>

                        if (.not. ma_push_get(mt_dbl,dimu*dime,
     &                             'c hole ue',l_uhole,k_uhole))
                        if (.not. ma_push_get(mt_dbl,dimv*dimf,
     &                             'c hole vf',l_vhole,k_vhole))
                        call ga_get(g_hole(espin),ulo,uhi,
     &                             v_sym(1,syme,espin),
     &                             v_sym(2,syme,espin),k_uhole,ld)
                        call ga_get(g_hole(fspin),vlo,vhi,
     &                             v_sym(1,symf,fspin),
     &                             v_sym(2,symf,fspin),k_fhole,ld)
                        call dgemm('n','n',dimv,dime,dimu,1.0d0,
     &                             dbl_mb(k_vulist),dimv,
     &                             dbl_mb(k_uhole),dimu,1.0d0,
     &                             dbl_mb(k_eflist),dimv)
                        call dgemm('t','n',dime,dimf,dimv,1.0d0,
     &                             dbl_mb(k_eflist),dimv,
     &                             dbl_mb(k_vhole),dimv,1.0d0,
     &                             dbl_mb(k_eflist),dime)
                        if (.not. ma_pop_stack(l_vhole)) call errquit()
                        if (.not. ma_pop_stack(l_uhole)) call errquit()
c                    I(me,nf)(symnf) -> R_temp(me,nf)(symnf) -> ga_copy_patch
                        call ga_put(g_r_array,
     &                             m_address,dime,dimm,
     &                             n_address,dimf,dimn,
     &                             dbl(k_eflist),dime,1.0d0)
                        if (.not. ma_pop_stack(l_eflist)) call errquit()
                        if (.not. ma_pop_stack(l_vulist)) call errquit()
c                    integral ga we want is lists(symmn,spin block)
c                    how is <ef|mn> stored symmetry wise
c                    STORED LIKE I(nv,mu) x for each irrep(mu) -> (um||vn) -> <uv||mn> 
c                    HENCE  LIKE I(nf,me) x for each irrep(me)
                     end do
                  end do
c           How is R stored on disk?  It is stored as the amplitudes
c           are in the triples. as R(m,e,symm,n,f,symf,symnf)
c           want to transform <uv|mn> to <ef|mn>=R(me,nf) -> R<fe|nm>
c           accumulate R to disk
               offset_me = file_t2_off(spinm,spine,spinn,spinf,symnf)+ov_off(1,symme,mspin,espin)
               offset_nf = file_t2_off(spinm,spine,spinn,spinf,symnf)+ov_off(1,symnf,nspin,fspin)
               call ga_copy_patch('n',g_r_array,1,dimme,1,dimnf,g_r,
     &                            offset_me,offset_me+dimme,offset_nf,offset_nf+dimnf)
               if (.not. ga_destroy(g_r_array)) call errquit()
            end do
         end do
      enddo
c
      end

