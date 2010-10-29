      subroutine uccsdt_cterm(urange,vrange,qO_handles,qV_handles)
c
c$Id$
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7), vrange(2,0:7)
      integer qO_handles(0:7,8), qV_handles(0:7,8)
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

c
c        generate T(iu,nf) aa,aa; ab,ab; ab,ba; ba,ab; ba,ba; bb,bb;
c        T stored as T(ia,jb,irrep) = T(i,a,symj,j,b,symb,symjb)
c        getT2 = uccsdt_ampfile_read_t2(file,spini,spina,spinj,spinb,
c                symjb,blo,bhi,g_t2,ocreate,distribution=block or column)
c     

c
c     Pure aa C terms (e=f=m=n=alpha)
c

C
C     GENERATION OF C-TERMS OF ALPHA-ALPHA SPIN-BLOCK. 
C     Z-INTERMEDIATES ARE USED TO GENERATE NECESSARY Q-INTERMEDIATES FOR E-TERMS
C
      do symnf = 0, nir-1
         symiu = symnf
         symme = symnf
         symjv = symiu
c
c     we need a subroutine to get the correct t-block and to do the transformation
c     collect z-terms with tjbnfmix and tjbnfpure
c     deallocate t-terms
c     generate q-terms
c     now generate tiumepure and tiumemix and combine with Z-terms into r-term
c     deallocate t-terms
c 
         dimnf = ov_len(symnf,1,1)
         dimiu_a = iu_len(symiu,1)
         dimiu_b = iu_len(symiu,2)
         dimjv_a = iso_len(symjv,1)
         dimjv_b = iso_len(symjv,2)
c        allocate TJBNFmix
c        get T(jb,nf,symnf) or T(jv,nf,symnf) j=b=beta n=f=alpha
         call get_T(g_t2_jbnf_mix,symnf,2,2,1,1,urange)
c
c     Z2(iu,nf) = <ij|ub>*T(jb,nf) i=u=alpha j=b=beta (integral list 2)
c     Z2(iu,nf) = I(jv,iu)^t * T(jv,nf)
c     Q2 from Z2
c     <if||un> += Z2 i=u=f=n=alpha (integral list 9 modified) = I(nf,iu)
c
         if (.not. ga_create(MT_DBL,dimnf,dimiu_a,'z2',0,0,g_z2) 
     &      call errquit('cterm: alloc z2 failed',0,0)
         call ga_dgemm('t','n',dimnf,dimiu_a,dimjv_b,1.0d0,
     &                 g_t2_jbnf_mix,list(symiu,2),0.0d0,g_z2)
         call make_qO(qO_handles(0,2),g_z2,symnf,1,1)
         call make_qV(qV_handles(0,2),g_z2,symnf,1,1)
         call oso_moindx2_add(g_z2,list(symnf,9),symnf,1,1,1,1)
c
c     Z3(iu,nf) = <ij||ub> * T(jb,nf)  i=u=j=b=beta  (integral list 4)
c     R(me,nf)  = T(iu,me)*[<fi|nu> + 0.5*Z3(iu,nf)]     (integral list 11)
c
         if (.not. ga_create(MT_DBL,dimnf,dimiu_b,'z3',0,0,g_z3) 
     &       call errquit('cterm: alloc z3 failed',0,0)
         call ga_dgemm('t','n',dimnf,dimiu_b,dimjv_b,1.0d0,
     &      g_t2_jbnf_mix,list(symiu,4),0.0d0,g_z3)
         if (.not. ga_destroy(g_t2_jbnf_mix)) call 
     &       errquit('cterm: dealloc g_t2_jbnf_mix failed',0)
         call ga_scale(g_z3,0.5d0)
         call oso_moindx2_add(g_z3,list(symnf,11),symnf,1,1,2,2)
c        MxM T*Z3 -> R(me,nf) all me, nf in block
c        get T(me,iu,symiu/symme)
         call get_T(g_t2_iume_mix,symiu,1,1,2,2,urange)
         call ga_dgemm('n','t',dimme,dimnf,dimiu_b,1.0d0,
     &      g_t2_iume_mix,g_z3,1.0d0,g_raa)
         if (.not. ga_destroy(g_t2_iume_mix)) call 
     &       errquit('cterm: dealloc g_t2_iume_mix failed',0)
         if (.not. ga_destroy(g_z3)) call
     &       errquit('cterm: dealloc g_z3 failed',0)

c     Here b could be v and T(jv,nf)
c     Z1(iu,nf) = <ij||ub>*T(jb,nf) i=j=u=b=n=f=alpha (integral list 1)
c     Q1
c     R(me,nf) = T(iu,me)*(<if||un>+Z2(iu,nf)+0.5*Z1(iu,nf))
c     <if||un> += Z2(iu,nf) (modified Z2)
c
         call get_T(g_t2_jbnf_pure,symnf,1,1,1,1,urange)
         if (.not. ga_create(MT_DBL,dimnf,dimiu_a,'z1',0,0,g_z1,
     &      call errquit('cterm: alloc z1 failed',0,0)
         call ga_dgemm('t','n',dimnf,dimiu_a,dimjb_a,1.0d0,
     &      g_t2_jbnf_pure,list(symiu,1),0.0d0,g_z1)
         call make_qO(qO_handles(0,1),g_z1,symnf,1,1)
         call make_qV(qV_handles(0,1),g_z1,symnf,1,1)
         call ga_add(0.5d0,g_z1,1.0d0,g_z2,g_z1)
         call ga_dgemm('n','n',dimnf,dimnf,dimiu_a,1.0d0,
     &      g_t2_iume_pure,g_z1,1.0d0,g_raa)
         if (.not. ga_destroy(g_t2_jbnf_pure)) call 
     &       errquit('cterm: dealloc g_t2_jbnf_pure failed',0)
         if (.not. ga_destroy(g_z2)) call 
     &       errquit('cterm: dealloc g_z2 failed',0)
         if (.not. ga_destroy(g_z1)) call 
     &       errquit('cterm: dealloc g_z1 failed',0)
      end do
C
C     GENERATION OF C-TERMS OF BETA-BETA SPIN-BLOCK PLUS
C     GENERATION OF C1- and C2-TERMS OF MIXED SPIN-BLOCK
C     Z-INTERMEDIATES ARE USED TO GENERATE NECESSARY Q-INTERMEDIATES FOR E-TERMS
C
      do symnf = 0, nir-1
         symiu = symnf
         symme = symnf
         symjv = symiu
         dimnf = ov_len(symnf,2,2)
         dimiu_a = iu_len(symiu,1)
         dimiu_b = iu_len(symiu,2)
         dimjv_a = iso_len(symjv,1)
         dimjv_b = iso_len(symjv,2)
c        allocate TIUMEpure
c        get T(iu,me) i=u=m=e=beta
         call get_T(g_t2_iume_pure,symme,2,2,2,2,urange)
c        allocate TIUMEmix
c        get T(iu,me) i=u=alpha, m=e=beta
         call get_T(g_t2_iume_mix,symme,1,1,2,2,urange)
c
c     Here b could be v and T(jv,nf)
c     Z5(iu,nf) = <ij|ub>*T(jb,nf) i=u=beta j=b=alpha (integral list 3)
         if (.not. ga_create(MT_DBL,dimnf,dimiu_b,'z5',0,0,
     &       g_z5) call errquit('cterm: alloc z5 failed',0,0)
         call ga_dgemm('t','n',dimnf,dimiu_a,dimjv_b,1.0d0,
     &                 g_t2_jbnf_mix,list(symiu,3),0.0d0,g_z5)
         call make_qO(qO_handles(0,5),g_z5,symnf,2,2)
         call make_qV(qV_handles(0,5),g_z5,symnf,2,2)
c     <if||un> += Z5(iu,nf) (modified list 12)
         call oso_moindx2_add(g_z5,list(symnf,12),symnf,2,2,2,2)
         if (.not. ga_destroy(list(symiu,12))) call
     &       errquit('cterm: dealloc list 12 failed',0)
c     Z6(iu,nf) = <ij||ub> * T(jb,nf)  i=u=j=b=alpha  (integral list 1)
         if (.not. ga_create(MT_DBL,dimnf,dimiu_a,'z6',0,0,g_z6)
     &       call errquit('cterm: alloc z6 failed',0,0)
         call ga_dgemm('t','n',dimnf,dimiu_a,dimjv_a,1.0d0,
     &      g_t2_jbnf_mix,list(symiu,1),0.0d0,g_z6)
         if (.not. ga_destroy(g_t2_jbnf_mix)) call
     &       errquit('cterm: dealloc g_t2_jbnf_mix failed',0)
         if (.not. ga_create(MT_DBL,dimnf,dimiu_a,'z7',0,0,g_z7)
     &       call errquit('cterm: alloc z7 failed',0,0)
c     copy z6 to z7 for use in mixed spin
         call ga_copy(g_z6,g_z7)
         call ga_scale(g_z6,0.5d0)
         call oso_moindx2_add(g_z6,list(symnf,10),symnf,2,2,1,1)
c        MxM T*Z3 -> R(me,nf) all me, nf in block
c        get T(me,iu,symiu/symme)
         call get_T(g_t2_iume_mix,symiu,1,1,2,2,urange)
         call ga_dgemm('n','t',dimme,dimnf,dimiu_a,1.0d0,
     &      g_t2_iume_mix,g_z6,1.0d0,g_rbb)
         if (.not. ga_destroy(g_z6)) call
     &       errquit('cterm: dealloc g_z6 failed',0)
         if (.not. ga_destroy(g_t2_iume_mix)) call
     &       errquit('cterm: dealloc g_t2_iume_mix failed',0)
c     Z4(iu,nf) = <ij||ub>*T(jb,nf) i=j=u=b=n=f=beta (integral list 4)
         call get_T(g_t2_jbnf_pure,symnf,1,1,1,1,urange)
         if (.not. ga_create(MT_DBL,dimnf,dimiu_a,'z4',0,0,g_z4,
     &       call errquit('cterm: alloc z4 failed',0,0)
         call ga_dgemm('t','n',dimnf,dimiu_a,dimjb_a,1.0d0,
     &      g_t2_jbnf_pure,list(symiu,4),0.0d0,g_z4)
         call make_qO(qO_handles(0,4),g_z4,symnf,2,2)
         call make_qV(qV_handles(0,4),g_z4,symnf,2,2)
         call ga_add(0.5d0,g_z4,1.0d0,g_z5,g_z5)
         call ga_dgemm('n','n',dimnf,dimnf,dimiu_a,1.0d0,
     &      g_t2_iume_pure,g_z5,1.0d0,g_rbb)
c     need the additional contribution of 1/2 Z4 for C2 term in mixed spin
         call ga_add(0.5d0,g_z4,1.0d0,g_z5,g_z5)
         if (.not. ga_destroy(g_t2_jbnf_pure)) call
     &       errquit('cterm: dealloc g_t2_jbnf_pure failed',0)
         if (.not. ga_destroy(g_z4)) call
     &       errquit('cterm: dealloc g_z4 failed',0)
         if (.not. ga_destroy(list(symiu,1))) call 
     &       errquit('cterm: dealloc list 1 failed',0)
         if (.not. ga_destroy(list(symiu,4))) call 
     &       errquit('cterm: dealloc list 4 failed',0)
C
C   WE WILL DO SOME MIXED TERMS AS WE ARE REUSING SOME BETA Z-BLOCKS
C
c 
c g_z5 (which is z4+z5+<if||an>) and g_z7 (which contains z6)
c are still active. We will use those in C2 and C1 mixed spin terms
c
c   
c   C2 += t(me,ia)([if|an> + Z4 + Z5] 
c   term in [] already stored in g_z5, just matmul with t
c
         call get_T(g_t2_meia_mix,symnf,1,1,2,2,urange)
         call ga_dgemm('n','n',dimnf,dimnf,dimiu_b,1.0d0,
     &      g_t2_iume_pure,g_z5,1.0d0,g_rab)
         if (.not. ga_destroy(g_t2_jbnf_pure)) call
     &       errquit('cterm: dealloc g_t2_jbnf_pure failed',0)
         if (.not. ga_destroy(g_z7)) call
     &       errquit('cterm: dealloc g_z1 failed',0)
c
c   Make Z7 with T is pure beta
c   Add to g_z6 and matmul with t
c   C1 += -t(ie,ma)*[<if|an> + Z6 + Z7] with T is pure alpha
c
         call get_T(g_t2_jbnf_pure,symnf,2,2,2,2,urange)
         call ga_dgemm('t','n',dimnf,dimiu_a,dimjv_b,1.0d0,
     &                 g_t2_jbnf_pure,list(symiu,2),1.0d0,g_z7)
         call oso_moindx2_add(g_z7,list(symnf,10),symnf,1,1,1,1)
         call get_T(g_t2_iema_pure,symnf,1,1,1,1,urange)
         call ga_dgemm('n','n',dimnf,dimnf,dimiu_a,1.0d0,
     &      g_t2_iume_pure,g_z7,1.0d0,g_rab)
         if (.not. ga_destroy(g_t2_iema_pure)) call
     &       errquit('cterm: dealloc g_t2_jbnf_pure failed',0)
         if (.not. ga_destroy(g_t2_jbnf_pure)) call
     &       errquit('cterm: dealloc g_t2_jbnf_pure failed',0)
         if (.not. ga_destroy(g_z7)) call
     &       errquit('cterm: dealloc g_z1 failed',0)
         if (.not. ga_destroy(list(symiu,2))) call 
     &       errquit('cterm: dealloc list 2 failed',0)
         if (.not. ga_destroy(list(symiu,10))) call 
     &       errquit('cterm: dealloc list 10 failed',0)
      end do
c
c
c     Mixed-spin terms (e=m=alpha n=f=beta)
c
C
C     GENERATION OF C-TERMS OF REMAINING MIXED ALPHA-BETA SPIN BLOCK
C
      do symnf = 0, nir-1
         symiu = symnf
         symme = symnf
         symjv = symiu
         dimnf = ov_len(symnf,2,2)
         dimnf = ov_len(symnf,2,2)
         dimiu_a = iu_len(symiu,1)
         dimiu_b = iu_len(symiu,2)
         dimjv_a = iso_len(symjv,1)
         dimjv_b = iso_len(symjv,2)
c
c        do blocks nf n=f=beta
c
c        allocate R(me,nf)
c
c     C4 : I DON'T THINK WE NEED THIS TERM, REARRANGEMENT REMOVED THIS ONE
c     Another piece of C5 using list 9
c     R(me,nf) += t(iu,nf)*[<ie||um> + Z2(iu,me)] 
c            I think we have the term inside [] in Z2, only multiply with T
c     .      e=m=i=u=alpha n=f=beta  (use modified list 9)
c
c        allocate TIUNFmix
c        get T(iu,nf) i=u=alpha, n=f=beta
         call get_T(g_t2_iunf_mix,symnf,1,1,2,2,urange)
         MxM I(iu,me)*T(iu,nf) -> R  (list 9)
         call ga_dgemm('t','n',dimiu,dimnf,dimjb,1.0d0,
     &      list(symiu,9),g_t2_iunf_mix,1.0d0,g_rab)
         free TIUNFmix
c
c     C5
c     R(me,nf) += t(iu,nf)*<ei|mu> 
c     .      e=m=alpha  n=f=i=u=beta (list 11)
c
c        allocate TIUNFpure
c        get T(iu,nf) i=u=n=f=beta
         call get_T(g_t2_iunf_pure,symnf,1,1,1,1,urange)
c        MxM I(iu,me)*T(iu,nf) -> R (list 11)
         call ga_dgemm('t','n',dimiu,dimnf,dimjb,1.0d0,
     &      list(symiu,11),g_t2_iunf_pure,1.0d0,g_rab)

c        accumulate R to disk
c        free R
         end do
c        free lists 9 and 11
         if (.not. ga_destroy(list(symiu,9))) call errquit()
         if (.not. ga_destroy(list(symiu,11))) call errquit()
      end do
c
c
      do symfn = 0, nir-1
         symme = symfn
         do symf = 0, nir-1
            symn = ieor(symfn,symf)
            do blocks of f
               allocate & zero R(me,nf)
               do symi = 0, nir-1
                  symif = ieor(symi,symf)
                  allocate T(mu,if) m=u=alpha i=f=beta
                  get T(mu,if) (read #sym times due to sym(fn) loop)
c     
c     C5
c     R(mf,ne) += - t(iu,mf)*<ei|un>
c                            I(en,iu)
c     .      e=m=u=alpha n=f=i=beta  (list 14)
c
c     Complexity here is that the OV pairs are mixed spin
c     whereas in most other uses each OV pair is pure spin
c     so we end up having to do an in-core transpose.
c
                  do symm = 0, nir-1
                     syme = ieor(symme,symm)
                     symu = syme
                     allocate R(mf,ne) ... dense 4-D array
                     allocate T(iu,mf) ... dense 4-D array
                     transpose T(mu,if) into T(iu,mf)
                     call ga_transpose(t_array)
                     MxM T(iu,mf)*I(iu,ne) -> R(mf,ne)
                     call ga_dgemm('t','n',dimiu,dimnf,dimjb,1.0d0,
     &                    list(symiu,14),g_t2_jbnf_pure,1.0d0,g_z4)
                     transpose accumulate R(mf,ne) into R(me,nf)
                  end do
                  free T(mu,if)
               end do
               accumulate R(me,nf) to disk
               free R(me,nf)
            end do
         end do
      end do
c
      destroy list(symnf,14)
c
c     Make Z8
c
c     Z8(iu,mf) = <ij|bu> * t(mb,jf) m=ib=alpha f=j=a=beta
c     .         = I(iu,jb) * T(jb,mf) (transposed list 3)
c
      allocate Z8
      do symjf = 0, nir-1
         symmb = symjf
         do symf = 0, nir-1
            symj = ieor(symjf,symf)
            do blocks of f
               allocate T(mb,jf)
               get T(mb,jf)
                  do symb = 0, nir-1
                  symm = ieor(symmb,symb)
                  


               end do
            end do
         end do
c        do sum_i C_hole(ua) Z8(iu,if) and store in q_8_O(a,f) -> qO_handles(8)
c        do sum_u C_hole(ua) Z8(iu,ma) and store in q_8_V(i,m) -> qV_handles(8)
         call make_qO(qO_handles(0,8),g_z8,symjf,0,1)
         call make_qV(qV_handles(0,8),g_z8,symjf,0,1)
         destroy list(symnf,3)
      end do
c
c   C3 += -t(ie,na)*[(if|ma) - Z8]
c

         destroy list(symnf,13)
      end
 
      subroutine get_T(g_t2,symjb,spini,spina,spinj,spinb,urange)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer spini,spina,spinj,spinb,symjb,g_t2
      integer blo,bhi,urange(2,0:7)
      integer t2_file
c
c        getT2 = uccsdt_ampfile_read_t2(file,spini,spina,spinj,spinb,
c                symjb,blo,bhi,g_t2,ocreate,distribution=block or column)
c
      blo = v_sym(1,0,spinb)
      bhi = v_sym(2,nir-1,spinb)
      getT2 = uccsdt_ampfile_read_t2(t2_file,spini,spina,spinj,spinb,
     &        symjb,blo,bhi,g_temp,.true.,'block')
c
c     allocate g_t2 with dimension for urange
c     do dgemm with inverse of c_hole to get u back
c
      end

      subroutine make_qO(qO_handles,g_z,symif,spini,spinf,urange)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer qO_handles(0:7),urange(1:2,0:7),g_z,symif,spini
      integer spinf,symf,symu,symi,len_i,iindx,jindx,u,f
      double precision fu_val, ifiu_val
      logical i_match,j_match
c
c     WE NEED TO ACCUMULATE THE Q's DURING THE BLOCKING LOOP IN acefterms.F
c     HENCE ADD CONTRIBUTIONS TO qO_handles
c     spini=spinj and spinf=spina
c     do sum_i C_hole(ua) Z2(iu,if) and store in q_2_O(a,f) -> qO_handles(2)
c     
      call ga_distribution(g_z,ga_nodeid(),ilo,ihi,jlo,jhi)
      jindx = 0
      do symf = 0, nir-1
         symu = symf
         symi = ieor(symif,symf)
         len_i = no(symi,spini)
         len_u = urange(2,symu)-urange(1,symu)+1
         len_f = v_sym(2,symf,spinf)-v_sym(1,symf,spinf)+1
         if (.not. ma_push_get(mt_dbl,len_f,'uf block',l_uf,
     &         k_uf)) call errquit('cterm: uf block alloc failed',0,0)
         do u = urange(1,symu), urange(2,symu)
            call ga_get(c_hole(spinf),urange(1,symu),urange(2,symu),
     &               v_sym(1,symf,spinf),v_sym(2,symf,spinf),
     &               dbl_mb(k_uf),len_u)
            iindx = oso_u_off(u,symif,spini)
            do f = v_sym(1,symf,spinf), v_sym(2,symf,spinf)
               i_match = ilo.le.(iindx+len_i).and.ihi.gt.(iindx)
               j_match = jlo.le.(jindx+len_i).and.jhi.gt.(jindx)
               if (i_match.and.j_match) then
                  do i = max(il0,iindx), min(ihi,iindx+len_i)
                     do j = max(jl0,jindx), min(jhi,jindx+len_i)
                        call ga_get(g_z,i,i,j,j,ifiu_val)
                        fu_val = fu_val + ifiu_val
                     end do
                  end do
                  call ga_acc(qO_handles(symf),alo,ahi,f,f,
      &                       dbl_mb(k_uf),len_f,fu_val)
               endif
               iindx = iindx + len_i
            end do
            jindx = jindx + len_i
         end do
      end do 
c 
      return
      end

      subroutine make_qV(qV_handles,g_z,symif,spini,spinf,urange)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer qV_handles(0:7),urange(1:2,0:7),g_z,symif,spini
      integer spinf,symf,symu,symi,len_i,iindx,jindx,u,f
      double precision fu_val, ifiu_val
      logical i_match,j_match
c
c     WE NEED TO ACCUMULATE THE Q's DURING THE BLOCKING LOOP IN acefterms.F
c     HENCE ADD CONTRIBUTIONS TO qV_handles
c     spini=spinj and spinf=spina
c     do sum_u C_hole(ua) Z2(ma,iu) and store in q_2_V(i,m) -> qV_handles(2)
c
      call ga_distribution(g_z,ga_nodeid(),ilo,ihi,jlo,jhi)
      jindx = 0
      do symf = 0, nir-1
         symu = symf
         symi = ieor(symif,symf)
         len_i = no(symi,spini)
         len_u = urange(2,symu)-urange(1,symu)+1
         len_f = v_sym(2,symf,spinf)-v_sym(1,symf,spinf)+1
         if (.not. ma_push_get(mt_dbl,len_i*len_i,'im block',l_im,
     &         k_im)) call errquit('cterm: im block alloc failed',0,0)
         if (.not. ma_push_get(mt_dbl,len_f,'uf block',l_uf,
     &         k_uf)) call errquit('cterm: uf block alloc failed',0,0)
         do u = urange(1,symu), urange(2,symu)
            call ga_get(c_hole(spinf),urange(1,symu),urange(2,symu),
     &               v_sym(1,symf,spinf),v_sym(2,symf,spinf),
     &               dbl_mb(k_uf),len_u)
            iindx = oso_u_off(u,symif,spini)
            do f = v_sym(1,symf,spinf), v_sym(2,symf,spinf)
               i_match = ilo.le.(iindx+len_i).and.ihi.gt.(iindx)
               j_match = jlo.le.(jindx+len_i).and.jhi.gt.(jindx)
               if (i_match.and.j_match) then 
                 alo = max(il0,iindx)
                 ahi = min(ihi,iindx+len_i)
                 mlo = max(jl0,jindx)
                 mhi = min(jhi,jindx+len_i)
                 call ga_get(g_z,alo,ahi,mlo,mhi,dbl_mb(k_im),ahi-alo+1)
                 call ga_acc(qV_handles(symf),alo,ahi,mlo,mhi,
      &                   dbl_mb(k_im),ahi-alo+1,dbl_mb(k_uf+fpointer)
               endif
               iindx = iindx + len_i
            end do
            jindx = jindx + len_i
         end do
      end do 
c  
      end


      subroutine oso_moindx2_add(g_a,g_b,sym,sa,sb,sc,sd)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer g_a,g_b,sym,sa,sb,sc,sd
      integer iu,symf,symn,dimf,dimn,dimv,dimiu
      integer l_nv,k_nv,l_nf,k_nf,l_vf,k_vf
      integer nvlo,nvhi,vlo,vhi,flo,fhi,nflo,nfhi
c
c     Do the transformation (nv,iu) -> (nf,iu) 
c     Transformation is for C_part
c     Hence nv -> nf or v -> f 
c     Loop over iu on the nodes and do a complete nv block at a time
c
      dimiu = iu_len(symiu,sc)
      do iu = ga_nodeid()+1, dimiu, ga_nnodes()
         do symf = 0, nir-1
            symn = ieor(sym,symf)
            dimf = nv_sym(symf,sb)
            dimn = no_sym(symn,sa)
            dimv = bf_per_ir(symf)
            if (.not. ma_push_get(mt_dbl,dimn*dimv,'nv block',l_nv,
     &         k_nv)) call errquit('cterm: nv block alloc failed',0,0)
            if (.not. ma_push_get(mt_dbl,dimn*dimf,'nf block',l_nf,
     &         k_nf)) call errquit('cterm: nf block alloc failed',0,0)
            if (.not. ma_push_get(mt_dbl,dimv*dimf,'vf block',l_vf,
     &         k_vf)) call errquit('cterm: vf block alloc failed',0,0)
            nvlo = oso_off(bf_per_ir_cum(symf)+1,sym,1) + 1
            nvhi = vlo + dimn*dimv
            vlo = bf_per_ir_cum(symf)+1
            vhi = vlo + dimv
            flo = v_sym(1,symf,sb)
            fhi = v_sym(2,symf,sb)
            call ga_get(g_b,vlo,vhi,iu,iu,dbl_mb(k_nv),dimn)
            call ga_get(g_part(sb),vlo,vhi,flo,fhi,dbl_mb(k_vf),dimv)
            call dgemm('n','n',dimn,dimf,dimv,1.0d0,dbl_mb(k_nv),
     &                 dimn,dbl_mb(k_vf),dimv,1.0d0,dbl_mb(k_nf),dimn)
            nflo = ov_off(v_sym(1,symf,1),0:7,sa,sb) + 1
            nfhi = flo + dimn*dimf
            call ga_acc(g_a,nflo,nfhi,iu,iu,dbl_mb(k_nf),dimn*dimf,
     &                  1.0d0)
            if (.not. ma_pop_stack(l_vf)) call
     &          errquit('oso_moindx2_add: vf block dealloc failed',0)
            if (.not. ma_pop_stack(l_nf)) call
     &          errquit('oso_moindx2_add: nf block dealloc failed',0)
            if (.not. ma_pop_stack(l_nv)) call
     &          errquit('oso_moindx2_add: nv block dealloc failed',0)
         end do
      end do
c
      end
