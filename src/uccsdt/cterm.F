      subroutine uccsdt_cterm()
      implicit none

c

c     Determine the blocking for SOs u in 
c     natural order ... and corresponding blocked range

      do blocks of SOs u
c
         generate 14 classes of partially transformed integrals
         that must be stored symmetry blocked.
c
c     i,j = particles
c     m,n = holes
c     a,b = holes
c
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
c     Lists 5,6,7,8 stored as <uv|mn> distributed over mn.
c
c     5. (mu|nf) m=alpha n=alpha       antisymmetrized
c     6. (mu|nf) m=alpha n=beta
c     7. (mu|nf) m=beta  n=alpha ... NOT NEEDED?
c     8. (mu|nf) m=beta  n=beta        antisymmetrized
c
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

c
c     Have lists 5-8 <uv|mn> distributed over mn ... transform
c     uv to ef, and accumlate <ef|mn> to disk.
c     <ef||mn> -> R(me,nf)
c
         do spin cases
            do sym(nf)
               -> sym(em)
               do sym(f)
                  -> sym(n) -> sym(m)
c     want to transform <uv|mn> to <ef|mn>=R
c
c     Assume for now that can hold all amplitudes for
c     given symmetry subblock????????
c
                  do blocks of f (hold R(me,nf) in memory)
                     allocate and zero R
                     do nm in parallel nm for given sym(n) and sym(m)
                        transform uv -> ef
                        put into R
                     end do
c     How is R stored on disk?  It is stored as the amplitudes
c     are in the triples.
                     accumulate R to disk
                  end do
               end do
            end do
         end do
c
c     <ij||mn> -> accumulate to GA where the corresponding X is stored.
c
         do spin cases
            do mn in parallel
               -> sym(ij)
               transform uv -> ij
            end do
            add in contribution to <ij||mn>
         end do
c
         free lists 5,6,7,8
c
         generate T(iu,nf) aa,aa; ab,ab; ab,ba; ba,ab; ba,ba; bb,bb;

c     

c$$$         do blocks of fn
c$$$            make Ziu,fn
c$$$            add Iiu,fn
c$$$            *Tiu,em
c$$$            store
c$$$         end do
                  
c
c     Pure aa C terms (e=f=m=n=alpha)
c
         do sym(fn) = sym(me)
c
            allocate TIUMEpure
            get T(iu,me) i=u=m=e=alpha
c
            allocate TIUMEmix
            get T(iu,me) i=u=beta, m=e=alpha
c
            do blocks nf n=f=alpha
c
n               allocate TJBNFmix
               get T(jb,nf) j=b=beta n=f=alpha
c
c     Z2(iu,nf) = <ij|ub>*T(jb,nf) i=u=alpha j=b=beta (integral list 2)
c     Q2
c     <if||un> += Z2 i=u=f=n=alpha (integral list 9 modified)
c
               allocate Z2
               MxM I*T -> Z2(iu,nf) i=u=n=f=alpha (integral list 2)
               Q2
               add Z2 to integral list 9
               free Z2
c
c     Z3(iu,nf) = <ij||ub> * T(jb,nf)  i=u=j=b=beta  (integral list 4)
c     R(me,nf)  = T(iu,me)*[<fi|nu> + Z3(iu,nf)]     (integral list 11)
c
               allocate Z3
               MxM I*T -> Z3(iu,nf) i=u=beta n=f=alpha
               add I(iu,nf)=<fi|nu>=(iu|nf) (list 3) to Z3
               MxM T*Z3 -> R(me,nf) all me, nf in block
               free Z3
               free TJBNFmix
c
c     Z1(iu,nf) = <ij||ub>*T(jb,nf) i=j=u=b=n=f=alpha (integral list 1)
c     Q1
c     R(me,nf) = T(iu,me)*(<if||un>+Z2(iu,nf)+0.5*Z1(iu,nf))
c     <if||un> += Z2(iu,nf) (modified list 9)
c
               allocate TJBNFpure
               get T(jb,nf) j=b=n=f=alpha
               allocate Z1
               MxM I*T -> Z1(iu,nf) 
               free TJBNFpure
               Q1
               scale Z1 by 0.5
               add modified list 9 to Z1
               MxM T*Z1 -> R(me,ne) all me, nf in block pure spin
               free Z1
c
            end do
            free TIUMEmix
            free TIUMEpure
         end do
c
c     Pure bb C terms (e=f=m=n=beta)
c
         do sym(fn) = sym(me)
c
            allocate TIUMEpure
            get T(iu,me) i=u=m=e=beta
c
            allocate TIUMEmix
            get T(iu,me) i=u=beta, m=e=beta
c
            do blocks nf n=f=beta
c
c     Z5(iu,nf) = <ij|ub>*T(jb,nf) i=u=beta j=b=alpha (integral list 3)
c     Q5
c     <if||un> += Z5(iu,nf) (modified list 12)
c
c     Z6(iu,nf) = <ij||ub> * T(jb,nf)  i=u=j=b=alpha  (integral list 1)
c     <if|un>  += Z6(iu,nf) (modified list 10)
c     R(me,nf)  = T(iu,me)*[<fi|nu> + Z6(iu,nf)]     (integral list 10)
c
c     Z4(iu,nf) = <ij||ub>*T(jb,nf) i=j=u=b=n=f=beta (integral list 4)
c     Q4
c
c     <if||un> += (1/2)*Z4(iu,nf) (modified list 12 again)
c     R(me,nf) = T(iu,me)*(<if||un>+Z5(iu,nf)+0.5*Z4(iu,nf)) (list 12)
c     <if||un> += (1/2)*Z4(iu,nf) (modified list 12 again)
c
            end do
         end do
c
         free lists 1 and 4
c
c     Mixed-spin terms (e=m=alpha n=f=beta)
c
         do sym(fn) = sym(me)
            do blocks nf n=f=beta
c
               allocate R(me,nf)
c
c     C4
c     R(me,nf) += t(iu,nf)*[<ie||um> + Z2(iu,me)] 
c     .      e=m=i=u=alpha n=f=beta  (use modified list 9)
c
               allocate TIUNFmix
               get T(iu,nf) i=u=alpha, n=f=beta
               MxM I(iu,me)*T(iu,nf) -> R  (list 9)
               free TIUNFmix
c
c     C5
c     R(me,nf) += t(iu,nf)*<ei|mu> 
c     .      e=m=alpha  n=f=i=u=beta (list 11)
c
               allocate TIUNFpure
               get T(iu,nf) i=u=n=f=beta
               MxM I(iu,me)*T(iu,nf) -> R (list 11)
c
               accumulate R to disk
               free R
            end do
         end do
c
         free lists 9 and 11
c
         do sym(fn) = sym(me)
            do sym(f)
               -> sym(n)
               do blocks of f
                  allocate & zero R(me,nf)
                  do sym(i)
                     -> sym(if)
                     allocate T(mu,if) m=u=alpha i=f=beta
                     get T(mu,if) (read #sym times due to sym(fn) loop)
c     
c     C6
c     R(mf,ne) += - t(iu,mf)*<ei|un>
c     .      e=m=u=alpha n=f=i=beta  (list 14)
c
c     Complexity here is that the OV pairs are mixed spin
c     whereas in most other uses each OV pair is pure spin
c     so we end up having to do an in-core transpose.
c
                     do sym(m)
                        -> sym(u)
                        -> sym(e)
                        allocate R(mf,ne) ... dense 4-D array
                        allocate T(iu,mf) ... dense 4-D array
                        transpose T(mu,if) into T(iu,mf)
                        MxM T(iu,mf)*I(iu,ne) -> R(mf,ne)
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
         do sym(jf)
            do sym(f)
               -> sym(j)
               do blocks of f
                  allocate T(mb,jf)
                  get T(mb,jf)
                  do sym(b)
                     -> sym(m)
                     
                  end do
               end do
            end do
         end do
      end do
c
      end
