      logical function jvltest(rtdb)
*     
*     $Id: jvltest.F,v 1.2 1998-10-28 02:25:29 jhvan Exp $
*      calling routine for CCSD(T) trials
*
      implicit none
      integer rtdb
c...    rtdb is handle to database
c
      call triple_so
c
      jvltest = .true.
c
      end
      subroutine triple_so
c
c...  first stab at a triple code for CCSD(T)
c...  using spin-orbital integrals sitting in core 
c...  for now a 4-dimensional array of <ij||ab> 
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
c
      integer occ,virt,nbasis
      logical oexist
      integer l_rint,k_rint,l_t1,k_t1,l_t2,k_t2,l_orb,k_orb
c
c...  check existence of required files
c
      inquire(file='T2',exist=oexist)
      if (.not. oexist)
     $     call errquit('jvltest: no T2 ', 0)
      inquire(file='ASOINTS',exist=oexist)
      if (.not. oexist)
     $     call errquit('jvltest: no ASOINTS ', 0)
      inquire(file='EVALS',exist=oexist)
      if (.not. oexist)
     $     call errquit('jvltest: no EVALS ', 0)
c
c... get basis size from T2
c
      open(9,file='T2',form='formatted')
      read(9,*) occ,virt
      nbasis = occ + virt
      close(9)
c
c..     integrals
c     double precision rint(nbasis,nbasis,nbasis,nbasis)
      if (.not. ma_push_get(mt_dbl, nbasis**4,'rint',l_rint, k_rint))
     $    call errquit('rint',0)
c..    amplitudes 
c     double precision t1(occ,virt), t2(occ,occ,virt,virt)
      if (.not. ma_push_get(mt_dbl, occ**2*virt**2,'t2',l_t2, k_t2))
     $    call errquit('t2',0)
      if (.not. ma_push_get(mt_dbl, occ*virt,'t1',l_t1, k_t1))
     $    call errquit('t1',0)
c..    the orb  (orbital energies)
c     double precision orb(nbasis)
      if (.not. ma_push_get(mt_dbl, nbasis,'orb',l_orb, k_orb))
     $    call errquit('orb',0)
c
      call get_INTS(dbl_mb(k_rint),nbasis)
      call get_T1(dbl_mb(k_t1),occ,virt)
      call get_T2(dbl_mb(k_t2),occ,virt)
      call get_EVALS(dbl_mb(k_orb),nbasis)
c
      call calc_robert(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
     1                dbl_mb(k_orb),nbasis,occ,virt)
c
      call calc_pople(dbl_mb(k_rint),dbl_mb(k_t1),dbl_mb(k_t2),
     1                dbl_mb(k_orb),nbasis,occ,virt)
c
c...  free all core beyond rint
      if (.not. ma_chop_stack(l_rint)) call errquit(' ma chop?', 0)
c
c      could also use ma_pop_stack(l_occ)  which is like gmem_free
c      (but using handle instead of pointer)
c
      end
      subroutine get_EVALS(orb,nbasis)
c
c...  real orbital energies from EVALS
c
      implicit none
c
      integer nbasis,nbin,i
      double precision orb(nbasis)
c
      integer nbin,i
      logical oexist
c
      inquire(file='EVALS',exist=oexist)
      if (.not.oexist) call errquit(' no EVAL ',0)
      open(9,file='EVALS',form='formatted')
      read(9,'(i5)') nbin
      if (nbin.ne.nbasis) call errquit('nbasis wrong',1)
      do i=1,nbasis
         read(9,'(5x,f20.10)') orb(i)
      end do
      close(9)
c
      print *,' orbital energies read '
c     print *, orb
      return
      end
      subroutine get_T2(t2,occ,virt)
c
c...  real t2 amplitudes from T2
c
      implicit none
c
      integer occ,virt,i
      double precision t2(occ,occ,virt,virt)
c
      integer ocin,vin,i,j,a,b
      logical oexist
c
      inquire(file='T2',exist=oexist)
      if (.not.oexist) call errquit(' no T2 ',0)
      open(9,file='T2',form='formatted')
      read(9,*) ocin,vin
      if (ocin.ne.occ.or.vin.ne.virt) call errquit(' occ,virt wrong',0)
c
      call dfill(occ*occ*virt*virt,0.0d0,t2,1)
c
10    read(9,'(4i5,f20.10)',end=20) i,j,a,b,t2(i,j,a,b)
      go to 10
c
20    close(9)
      print *,' T2 read '
c     print *,((((t2(i,j,a,b),i=1,3),j=1,3),a=1,3),b=1,3)
c
      return
      end
      subroutine get_T1(t1,occ,virt)
c
c...  real t1 amplitudes from T1
c
      implicit none
c
      integer occ,virt,i
      double precision t1(occ,virt)
c
      integer ocin,vin,i,a
      logical oexist
c
      inquire(file='T1',exist=oexist)
      if (.not.oexist) then
         print *,' no T1 '
         call dfill(occ*virt,0.0d0,t1,1)
         return
      end if
      open(9,file='T1',form='formatted')
      read(9,'(i5,i3)') ocin,vin
      if (ocin.ne.occ.or.vin.ne.virt) call errquit(' occ,virt wrong',1)
c
      call dfill(occ*virt,0.0d0,t1,1)
c
10    read(9,'(2i5,f20.10)',end=20) i,a,t1(i,a)
      go to 10
c
20    close(9)
      print *,' T1 read '
c
      return
      end
      subroutine get_INTS(rint,nbasis)
c
c...  real integrals from ASOINTS
c
      implicit none
c
      integer nbasis
      double precision rint(nbasis,nbasis,nbasis,nbasis)
c
      integer in,i,j,k,l
      logical oexist
c
      inquire(file='ASOINTS',exist=oexist)
      if (.not.oexist) call errquit(' no INTS ',0)
      open(9,file='ASOINTS',form='formatted')
      read(9,'(i5)') in
      if (in.ne.nbasis) call errquit(' nbasis on INTS wrong',0)
c
      call dfill(nbasis**4,0.0d0,rint,1)
c
10    read(9,'(4i5,f20.10)',end=20) i,j,k,l,rint(i,j,k,l)
      go to 10
c
20    close(9)
      print *,' integrals read '
c     print *,((((rint(i,j,k,l),i=1,3),j=1,3),k=1,3),l=1,3)
c
      return
      end
      subroutine calc_robert(rint,t1,t2,orb,nbasis,occ,virt)
c
c...   use Robert's formulae
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t1(occ,virt), t2(occ,occ,virt,virt)
      double precision orb(nbasis)
c..     intermediates
      double precision w1,w2,w3,w4,v,w1w1,w1w2,w1w3,w1w4
      double precision vw1,vw2,vw2,vw3,vw4
c
      double precision e1,e2,energy,delta
      integer i,j,k,a,b,c,e,m
c
      w1w1 = 0.0d0
      w1w2 = 0.0d0
      w1w3 = 0.0d0
      w1w4 = 0.0d0
      vw1 = 0.0d0
      vw2 = 0.0d0
      vw3 = 0.0d0
      vw4 = 0.0d0
c
      do i=1,occ
       do j=1,occ
        do k=1,occ
         do a=1,virt
          do b=1,virt
           do c=1,virt
            w1 = 0.0d0
            w2 = 0.0d0
            w3 = 0.0d0
            w4 = 0.0d0
            do e=1,virt
             w1 = w1 + t2(i,j,c,e)*rint(a+occ,b+occ,e+occ,k)
             w2 = w2 + t2(k,j,c,e)*rint(a+occ,b+occ,e+occ,i)
             w3 = w3 + t2(i,j,a,e)*rint(c+occ,b+occ,e+occ,k)
             w4 = w4 + t2(k,j,a,e)*rint(c+occ,b+occ,e+occ,i)
            end do
            do m=1,occ
             w1 = w1 + t2(k,m,a,b)*rint(c+occ,m,i,j)
             w2 = w2 + t2(i,m,a,b)*rint(c+occ,m,k,j)
             w3 = w3 + t2(k,m,c,b)*rint(a+occ,m,i,j)
             w4 = w4 + t2(i,m,c,b)*rint(a+occ,m,k,j)
            end do
            v = t1(k,c)*rint(i,j,a,b)
            delta = (orb(i)+orb(j)+orb(k)
     1            - orb(a+occ)-orb(b+occ)-orb(c+occ))
            w1w1 = w1w1 + w1*w1 / delta
            w1w2 = w1w2 + w1*w2 / delta
            w1w3 = w1w3 + w1*w3 / delta
            w1w4 = w1w4 + w1*w4 / delta
            vw1 = vw1 + v*w1 / delta
            vw2 = vw2 + v*w2 / delta
            vw3 = vw3 + v*w3 / delta
            vw4 = vw4 + v*w4 / delta
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  calculate energies
c
      e1 = (vw1 - 2.0d0*vw2 -2.0d0*vw3 + 4.0d0*vw4) /4.0d0
      e2 = (w1w1 - 2.0d0*w1w2 -2.0d0*w1w3 + 4.0d0*w1w4) /4.0d0
      energy = e1 + e2
c
      print *,' Robert : e1= ',e1,' e2= ',e2,' e= ',energy
c
      end
      subroutine calc_pople(rint,t1,t2,orb,nbasis,occ,virt)
c
c...   use Pople's formulae (who can argue with a Nobel-prize)
c...   BUT ... typed in from Laplace triple document (jvl)
c...   ( + formulae(QCI ...))
c
      implicit none
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t1(occ,virt), t2(occ,occ,virt,virt)
      double precision orb(nbasis)
c..     intermediates
      double precision u1,u2,u3,u4,u5,u6,u7,u8,u9,v1
c
      double precision e1,e2,energy,delta
      integer i,j,k,a,b,c,e,m
c
      e1 = 0.0d0
      e2 = 0.0d0
c
      do i=1,occ
       do j=1,occ
        do k=1,occ
         do a=1,virt
          do b=1,virt
           do c=1,virt
            u1 = 0.0d0
            u2 = 0.0d0
            u3 = 0.0d0
            u4 = 0.0d0
            u5 = 0.0d0
            u6 = 0.0d0
            u7 = 0.0d0
            u8 = 0.0d0
            u9 = 0.0d0
            do e=1,virt
             u1 = u1 + t2(i,j,a,e)*rint(b+occ,c+occ,e+occ,k)
             u2 = u2 + t2(i,j,b,e)*rint(c+occ,a+occ,e+occ,k)
             u3 = u3 + t2(i,j,c,e)*rint(a+occ,b+occ,e+occ,k)
             u4 = u4 + t2(k,i,a,e)*rint(b+occ,c+occ,e+occ,j)
             u5 = u5 + t2(k,i,b,e)*rint(c+occ,a+occ,e+occ,j)
             u6 = u6 + t2(k,i,c,e)*rint(a+occ,b+occ,e+occ,j)
             u7 = u7 + t2(j,k,a,e)*rint(b+occ,c+occ,e+occ,i)
             u8 = u8 + t2(j,k,b,e)*rint(c+occ,a+occ,e+occ,i)
             u9 = u9 + t2(j,k,c,e)*rint(a+occ,b+occ,e+occ,i)
            end do
            do m=1,occ
             u1 = u1 + t2(i,m,a,b)*rint(c+occ,m,j,k)
             u2 = u2 + t2(i,m,b,c)*rint(a+occ,m,j,k)
             u3 = u3 + t2(i,m,c,a)*rint(b+occ,m,j,k)
             u4 = u4 + t2(j,m,a,b)*rint(c+occ,m,k,i)
             u5 = u5 + t2(j,m,b,c)*rint(a+occ,m,k,i)
             u6 = u6 + t2(j,m,c,a)*rint(b+occ,m,k,i)
             u7 = u7 + t2(k,m,a,b)*rint(c+occ,m,i,j)
             u8 = u8 + t2(k,m,b,c)*rint(a+occ,m,i,j)
             u9 = u9 + t2(k,m,c,a)*rint(b+occ,m,i,j)
            end do
            v1 = t1(i,a)*rint(j,k,b+occ,c+occ) 
     1         + t1(i,b)*rint(j,k,c+occ,a+occ) 
     2         + t1(i,c)*rint(j,k,a+occ,b+occ) 
     3         + t1(j,a)*rint(k,i,b+occ,c+occ) 
     4         + t1(j,b)*rint(k,i,c+occ,a+occ)
     5         + t1(j,c)*rint(k,i,a+occ,b+occ) 
     6         + t1(k,a)*rint(i,j,b+occ,c+occ) 
     7         + t1(k,b)*rint(i,j,c+occ,a+occ)
     8         + t1(k,c)*rint(i,j,a+occ,b+occ) 
            delta = (orb(i)+orb(j)+orb(k)
     1            - orb(a+occ)-orb(b+occ)-orb(c+occ))
            e1 = e1+ v1*(u1+u2+u3+u4+u5+u6+u7+u8+u9) / delta
            e2 = e2+ (u1+u2+u3+u4+u5+u6+u7+u8+u9)**2 / delta
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  calculate energies
c
      e1 = e1 / 36.0d0
      e2 = e2 / 36.0d0
      energy = e1 + e2
c
      print *,' Pople  : e1= ',e1,' e2= ',e2,' e= ',energy
c
      end
#ifdef cullen
      subroutine calc_cullen(rint,t1,t2,orb,nbasis,occ,virt
     1                       wab,)
c
c...   use Cullen+Zerner's approach
c...   starting from "Robert's" formula
c...   approach 1 generate all n**4 intermediates
c
      implicit none
c
      dimension wab(occ,occ,virt,virt),wij(occ,occ,virt,virt)
      dimension wtab(occ,occ,virt,virt),wtij(occ,occ,virt,virt)
c
      integer nbasis,occ,virt
      double precision rint(nbasis,nbasis,nbasis,nbasis)
      double precision t1(occ,virt), t2(occ,occ,virt,virt)
      double precision orb(nbasis)
c..     intermediates
      double precision u1,u2,u3,u4,u5,u6,u7,u8,u9,v1
c
      double precision e1,e2,energy,delta
      integer i,j,k,a,b,c,e,m
c
c..     wab(occ,occ,virt,virt)=sum(ab)(ab||ei)*(ab||fj)
c..     wtab(occ,occ,virt,virt)=sum(ab)t(ij,ab)*(ab||fj)
c
      do e=1,virt
       do f=1,virt
        do i=1,occ
         do j=1,occ
          wab(i,j,e,f) = 0.0d0
          do a=1,virt
           do b=1,virt
            wab(i,j,e,f) = rint(a+occ,b+occ,e+occ,i)
     1                   * rint(a+occ,b+occ,f+occ,j)
           end do
          end do
         end do
        end do
       end do
      end do
c..     wij(occ,occ,virt,virt)=sum(ij) (an||ij)*(bm||ij)
      do 


      do e=1,virt
       do f=1,virt
        do k=1,occ
         wekfk = 0.0d0
         do a=1,nvirt
          do b=1,nvirt
           wekfk = wekfk + rint(a,b,e,k)*rint(a,b,f,k)
          end do
         end do
         do i=1,occ
          do j=1,occ
           do c=1,virt
            tijcef = tijcef + t2(i,j,c,e)*t2(i,j,c,f) nopt rihgt
       end do
      end do
      do i=1,occ
       do j=1,occ
        do k=1,occ
         do a=1,virt
          do b=1,virt
           do c=1,virt
            u1 = 0.0d0
            u2 = 0.0d0
            u3 = 0.0d0
            u4 = 0.0d0
            u5 = 0.0d0
            u6 = 0.0d0
            u7 = 0.0d0
            u8 = 0.0d0
            u9 = 0.0d0
            do e=1,virt
             u1 = u1 + t2(i,j,a,e)*rint(b,c,e,k)
             u2 = u2 + t2(i,j,b,e)*rint(c,a,e,k)
             u3 = u3 + t2(i,j,c,e)*rint(a,b,e,k)
             u4 = u4 + t2(k,i,a,e)*rint(b,c,e,j)
             u5 = u5 + t2(k,i,b,e)*rint(c,a,e,j)
             u6 = u6 + t2(k,i,c,e)*rint(a,b,e,j)
             u7 = u7 + t2(j,k,a,e)*rint(b,c,e,i)
             u8 = u8 + t2(j,k,b,e)*rint(c,a,e,i)
             u9 = u9 + t2(j,k,c,e)*rint(a,b,e,i)
            end do
            do m=1,occ
             u1 = u1 + t2(i,m,a,b)*rint(c,m,j,k)
             u2 = u2 + t2(i,m,b,c)*rint(a,m,j,k)
             u3 = u3 + t2(i,m,c,a)*rint(b,m,j,k)
             u4 = u4 + t2(j,m,a,b)*rint(c,m,k,i)
             u5 = u5 + t2(j,m,b,c)*rint(a,m,k,i)
             u6 = u6 + t2(j,m,c,a)*rint(b,m,k,i)
             u7 = u7 + t2(k,m,a,b)*rint(c,m,i,j)
             u8 = u8 + t2(k,m,b,c)*rint(a,m,i,j)
             u9 = u9 + t2(k,m,c,a)*rint(b,m,i,j)
            end do
            v1 = t1(i,a)*rint(j,k,b,c) + t1(i,b)*rint(j,k,c,a) 
     1         + t1(i,c)*rint(j,k,a,b) 
     2         + t1(j,a)*rint(k,i,b,c) + t1(j,b)*rint(k,i,c,a)
     3         + t1(j,c)*rint(k,i,a,b) 
     2         + t1(k,a)*rint(i,j,b,c) + t1(k,b)*rint(i,j,c,a)
     5         + t1(k,c)*rint(i,j,a,b) 
            delta = (orb(i)+orb(j)+orb(k)-orb(a)-orb(b)-orb(c))
            e1 = e1+ v1*(u1+u2+u3+u4+u5+u6+u7+u8+u9) / delta
            e2 = e2+ (u1*u1+u2*u2+u3*u3+u4*u4+u5*u5+u6*u6+u7*u7
     1           +u8*u8+u9*u9) / delta
           end do
          end do
         end do
        end do
       end do
      end do
c
c...  calculate energies
c
      e1 = e1 / 36.0d0
      e2 = e2 / 36.0d0
      energy = e1 + e2
c
      print *,' Pople  : e1= ',e1,' e2= ',e2,' e= ',energy
c
      end
#endif
