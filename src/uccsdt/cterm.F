      subroutine uccsdt_cterm(urange,vrange,qO_handles,qV_handles)
c
c$Id: cterm.F,v 1.3 2002-08-28 22:16:19 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer urange(2,0:7), vrange(2,0:7)
      integer qO_handles(8,0:7), qV_handles(8,0:7)
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
         dimiu_a = iu_len(symiu,0)
         dimiu_b = iu_len(symiu,1)
         dimjv_a = iso_len(symjv,0)
         dimjv_b = iso_len(symjv,1)
c        allocate TJBNFmix
c        if (.not. ga_create(MT_DBL,dimjv_b,dimnf,'Tjvnf mix',0,0,
c    &       g_t2_jbnf_mix) call errquit()
c        get T(jb,nf,symnf) or T(jv,nf,symnf) j=b=beta n=f=alpha
c        subroutine get_T(symcd,spina,spinb,spinc,spind) does get and transforms one index back to SO basis
         call get_T(g_t2_jbnf_mix,symnf,1,1,0,0,urange)
c
c     Here b could be v and T(jv,nf)
c     Z2(iu,nf) = <ij|ub>*T(jb,nf) i=u=alpha j=b=beta (integral list 2)
c     Z2(iu,nf) = I(jv,iu)^t * T(jv,nf)
c     Q2 from Z2
c     <if||un> += Z2 i=u=f=n=alpha (integral list 9 modified) = I(nf,iu)
c
         if (.not. ga_create(MT_DBL,dimiu_a,dimnf,'z2',0,0,g_z2) 
     &      call errquit()
         call ga_dgemm('t','n',dimiu_a,dimnf,dimjv_b,1.0d0,
     &                 list(symiu,2),g_t2_jbnf_mix,1.0d0,g_z2)
c        do sum_i C_hole(ua) Z2(iu,if) and store in q_2_O(a,f) -> qO_handles(2)
c           do sum over i first to get u,f and then do dgemm with c_hole
c        do sum_u C_hole(ua) Z2(iu,ma) and store in q_2_V(i,m) -> qV_handles(2)
c           need to do dgemm with c_hole first and then do sum
c            do both of them local on a node with get's puts and local sums
         call make_qO(qO_handles(2,0),g_z2,symnf,0,0)
         call make_qV(qV_handles(2,0),g_z2,symnf,0,0)
c        add Z2 to integral list 9, But list 9 is (nv,iu) !!!!!
c        transform v index in list 9 and transpose matrix
c        add Z2 to new transformed and transposed list 9
         call ga_add(1.0d0,g_z2,1.0,list(symiu,9),list(symiu,9)) 
c        no, transform v index in list(symiu,9) and add to Z2,
c        delete list 9 and use Z2 in aaC1 and abC4
         if (.not. ga_destroy(list(symiu,9))) call errquit()
c
c     Here b could be v and T(jv,nf)
c     Z3(iu,nf) = <ij||ub> * T(jb,nf)  i=u=j=b=beta  (integral list 4)
c     R(me,nf)  = T(iu,me)*[<fi|nu> + Z3(iu,nf)]     (integral list 11)
c
c     allocate Z3
         if (.not. ga_create(MT_DBL,dimiu_b,dimnf,'z3',0,0,g_z3) 
     &       call errquit()
c        MxM I*T -> Z3(iu,nf) i=u=beta n=f=alpha
         call ga_dgemm('t','n',dimiu_b,dimnf,dimjv_b,1.0d0,
     &      list(symiu,4),g_t2_jbnf_mix,1.0d0,g_z3)
c        free t2_jbnf_mix
         if (.not. ga_destroy(g_t2_jbnf_mix)) call errquit()
c        add I(nf,iu)=<fi|nu>=(iu|nf) (list 11) to Z3 after transforming third index v->f
c        transform v index in list 11 and transpose matrix and add to Z3
c        because of strided data get all iu for one nf per procesor
c        do a local dgemm to transform and add it to g_z3
         call ga_add(1.0d0,g_z3,1.0d0,list(symiu,11),g_z3) 
c        MxM T*Z3 -> R(me,nf) all me, nf in block
c        first create t2_iume_mix
c        if (.not. ga_create(MT_DBL,dimiu_b,dimme,'Tiume mix',0,0,
c    &       g_t2_iume_mix) call errquit()
c        get T(iu,me,symme) 
         call get_T(g_t2_iume_mix,symme,1,1,0,0,urange)
         call ga_dgemm('t','n',dimme,dimnf,dimiu_b,1.0d0,
     &      g_t2_iume_mix,g_z3,1.0d0,g_raa)
c        free t2_iume_mix
         if (.not. ga_destroy(g_t2_iume_mix)) call errquit()
c        free Z3
         if (.not. ga_destroy(g_z3)) call errquit()

c     Here b could be v and T(jv,nf)
c     Z1(iu,nf) = <ij||ub>*T(jb,nf) i=j=u=b=n=f=alpha (integral list 1)
c     Q1
c     R(me,nf) = T(iu,me)*(<if||un>+Z2(iu,nf)+0.5*Z1(iu,nf))
c     <if||un> += Z2(iu,nf) (modified list 9)
c
c     allocate TJBNFpure
c     get T(jb,nf) j=b=n=f=alpha
         call get_T(g_t2_jbnf_pure,symnf,0,0,0,0,urange)
c     allocate Z1
         if (.not. ga_create(MT_DBL,dimiu,dimnf,'z1',blockiu,
     &      blocknf,g_z1) call errquit()
c        MxM I*T -> Z1(iu,nf) 
         call ga_dgemm('t','n',dimiu,dimnf,dimjb,1.0d0,
     &      list(symiu,1),g_t2_jbnf_pure,1.0d0,g_z1)
c        free TJBNFpure
         if (.not. ga_destroy(g_t2_jbnf_pure)) call errquit()
c        do sum_i C_hole(ua) Z1(iu,if) and store in q_1_O(a,f) -> qO_handles(1)
c           do sum over i first to get u,f and then do dgemm with c_hole
c        do sum_u C_hole(ua) Z1(iu,ma) and store in q_1_V(i,m) -> qV_handles(1)
c           need to do dgemm with c_hole first and then do sum
c            do both of them local on a node with get's puts and local sums
         call make_qO(qO_handles(1,0),g_z1,symnf,0,0)
         call make_qV(qV_handles(1,0),g_z1,symnf,0,0)
c        scale Z1 by 0.5
         call ga_scale(g_z1,0.5d0)
c        add Z2 to Z1
         call ga_add(1.0d0,g_z1,1.0d0,g_z2,g_z1)
c        MxM T*Z1 -> R(me,nf) all me, nf in block pure spin
         call ga_dgemm('t','n',dimme,dimnf,dimiu,1.0d0,
     &      g_t2_iume_pure,g_z1,1.0d0,g_raa)
c        free Z1
         if (.not. ga_destroy(g_z1)) call errquit()
c        free t2_iume_pure
         if (.not. ga_destroy(g_t2_iume_pure)) call errquit()
      end do
c
c     Pure bb C terms (e=f=m=n=beta)
c
C
C     GENERATION OF C-TERMS OF BETA-BETA SPIN-BLOCK. 
C     Z-INTERMEDIATES ARE USED TO GENERATE NECESSARY Q-INTERMEDIATES FOR E-TERMS
C
      do symnf = 0, nir-1
         symme = symnf
         symiu = symnf

c        allocate TIUMEpure
c        get T(iu,me) i=u=m=e=beta
         call get_T(g_t2_iume_pure,symme,1,1,1,1,urange)
c        allocate TIUMEmix
c        get T(iu,me) i=u=alpha, m=e=beta
         call get_T(g_t2_iume_mix,symme,0,0,1,1,urange)
c
c        do blocks nf n=f=beta

c     Here b could be v and T(jv,nf)
c     Z5(iu,nf) = <ij|ub>*T(jb,nf) i=u=beta j=b=alpha (integral list 3)
            if (.not. ga_create(MT_DBL,dimiu,dimnf,'z5',blockiu,
     &         blocknf,g_z5) call errquit()
            call ga_dgemm('t','n',dimiu,dimnf,dimjb,1.0d0,
     &         list(symiu,3),g_t2_jbnf_mix,1.0d0,g_z5)
c        do sum_i C_hole(ua) Z1(iu,if) and store in q_1_O(a,f) -> qO_handles(5)
c        do sum_u C_hole(ua) Z1(iu,ma) and store in q_1_V(i,m) -> qV_handles(5)
         call make_qO(qO_handles(5,0),g_z1,symnf,1,1)
         call make_qV(qV_handles(5,0),g_z1,symnf,1,1)
c
c     <if||un> += Z5(iu,nf) (modified list 12)
            call ga_add(1.0d0,g_z5,1.0,list(symiu,12),list(symiu,12)) 
c
c                 Here b could be v and T(jv,nf)
c     Z6(iu,nf) = <ij||ub> * T(jb,nf)  i=u=j=b=alpha  (integral list 1)
            if (.not. ga_create(MT_DBL,dimiu,dimnf,'z6',blockiu,
     &         blocknf,g_z6) call errquit()
            call ga_dgemm('t','n',dimiu,dimnf,dimjb,1.0d0,
     &         list(symiu,1),g_t2_jbnf_mix,1.0d0,g_z6)
c     <if|un>  += Z6(iu,nf) (modified list 10)
            call ga_add(1.0d0,g_z6,1.0,list(symiu,10),list(symiu,10)) 
c     R(me,nf)  = T(iu,me)*[<fi|nu> + Z6(iu,nf)]     (integral list 10)
            call ga_dgemm('t','n',dimme,dimnf,dimiu,1.0d0,
     &         g_t2_iume_mix,list(symiu,10),1.0d0,g_rbb)
c
c                 Here b could be v and T(jv,nf)
c     Z4(iu,nf) = <ij||ub>*T(jb,nf) i=j=u=b=n=f=beta (integral list 4)
            if (.not. ga_create(MT_DBL,dimiu,dimnf,'z4',blockiu,
     &         blocknf,g_z4) call errquit()
            call ga_dgemm('t','n',dimiu,dimnf,dimjb,1.0d0,
     &         list(symiu,4),g_t2_jbnf_pure,1.0d0,g_z4)
c        do sum_i C_hole(ua) Z4(iu,if) and store in q_4_O(a,f) -> qO_handles(4)
c        do sum_u C_hole(ua) Z4(iu,ma) and store in q_4_V(i,m) -> qV_handles(4)
         call make_qO(qO_handles(4,0),g_z1,symnf,1,1)
         call make_qV(qV_handles(4,0),g_z1,symnf,1,1)
c
c     <if||un> += (1/2)*Z4(iu,nf) (modified list 12 again)
            call ga_add(0.5d0,g_z4,1.0,list(symiu,12),list(symiu,12)) 
c
c     R(me,nf) = T(iu,me)*(<if||un>+Z5(iu,nf)+0.5*Z4(iu,nf)) (list 12)
            call ga_dgemm('t','n',dimme,dimnf,dimiu,1.0d0,
     &         g_t2_iume_pure,list(symiu,4),1.0d0,g_rbb)
c     <if||un> += (1/2)*Z4(iu,nf) (modified list 12 again) to be used in mixed spin
               call ga_add(0.5d0,g_z4,1.0,list(symiu,12),list(symiu,12)) 
c
c
c        end do
c        free lists 1 and 4
         if (.not. ga_destroy(list(symiu,1))) call errquit()
         if (.not. ga_destroy(list(symiu,4))) call errquit()
      end do
c
c
c     Mixed-spin terms (e=m=alpha n=f=beta)
c
C
C     GENERATION OF C-TERMS OF MIXED ALPHA-BETA SPIN BLOCK
C
      do symnf = 0, nir-1
         symme = symnf
         symiu = symnf 
c        do blocks nf n=f=beta
c
         allocate R(me,nf)
c
c     C4
c     R(me,nf) += t(iu,nf)*[<ie||um> + Z2(iu,me)] 
c            I think we have the term inside [] in Z2, only multiply with T
c     .      e=m=i=u=alpha n=f=beta  (use modified list 9)
c
c        allocate TIUNFmix
c        get T(iu,nf) i=u=alpha, n=f=beta
         call get_T(g_t2_iunf_mix,symnf,0,0,1,1,urange)
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
      free list 14
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
         call make_qO(qO_handles(8,0),g_z1,symjf,0,1)
         call make_qV(qV_handles(8,0),g_z1,symjf,0,1)
      end do
c
c   C3 += -t(ie,na)*[(if|ma) - Z8]
c

c
c   Make Z7
c   C1 += -t(ie,ma)*[<if|an> + Z6 + Z7]
c


c   
c   C2 += t(me,ia)([if|an> + Z4 + Z5]
c

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

      subroutine make_qO(qO_handles,g_z,symif,spini,spinf)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer qO_handles(0:7)
c
c     spini=spinj and spinf=spina
c     do sum_i C_hole(ua) Z2(iu,if) and store in q_2_O(a,f) -> qO_handles(2)
c     do sum over i first to get u,f and then do dgemm with c_hole
c
c     symmetry wise i=i and hence u=f as iu=if, hence
c     do symf = 0, nir-1
c        symu = syma = symf
c        symi = ieor(symnf=symif,symf)
c        do loop over f
c            do loop over a
c               sum iho,ihi (symi) and transform with c_hole(u,a)^t * Z2(u,f)
c        

c 
      end

      subroutine make_qV(qV_handles,g_z,symif,spini,spinf)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer qV_handles(0:7)
c
c     spini=spinj and spinf=spina
c     do sum_u C_hole(ua) Z2(iu,ma) and store in q_2_V(i,m) -> qV_handles(2)
c     need to do dgemm with c_hole first and then do sum
c     do both of them local on a node with get's puts and local sums
c
c     symmetry wise i=i and hence u=f as iu=if, hence
c     do symf = 0, nir-1
c        symu = syma = symf
c        symi = ieor(symnf=symif,symf)
c        do loop over f
c            do loop over a
c               sum iho,ihi (symi) and transform with c_hole(u,a)^t * Z2(u,f)
c

c  
      end


