C> \file makehd.F
C> Make diagonal block of Hamiltonian
C>
C> \ingroup selci
C> @{
C>
C> \brief Make diagonal block of Hamiltonian
C>
C> Form the hamiltonian matrix within a spin block for fixed
C> orbital occupancy ... full square is generated in hii
C>
C> If odonly is true then only the actual diagonal of hii
C> is generated ... to optimize select
C>
      subroutine selci_makehd(hii,odonly,
     $     ns,nd,indbar,iocc,lists,listd,w1,w2,
     $              work1,work2,f,h,g,int12,int34,numf)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
      dimension w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1),
     $     indbar(*),iocc(*),work1(*),work2(*),hii(*),
     $     f(*),h(*),g(*),lists(*),listd(*),int12(*),int34(*)
      logical odonly
#include "stfunc.fh"
c
c     Form the hamiltonian matrix within a spin block for fixed
c     orbital occupancy ... full square is generated in hii
c
c     If odonly is true then only the actual diagonal of hii
c     is generated ... to optimize select
c
c     <Iu|H|Iv> = 1/2 { sum(i)[ni*(hii+fii)] + 
c                       sum(i docc,j socc)(ij|ji) } delta(u,v) +
c                 sum(i socc>j socc)(ij|ji)<u|EijEji|v>
c
c     fij = hij + sum(r)[nr * (ij|rr) - (ir|rj)]
c
      numf   = nf(ns)
      numfsq = numf*numf
      numf1  = numf + 1
      if (.not. odonly) call dfill(numfsq,0.0d0,hii,1)
c
c     do the true diagonal bit first
c
      e = 0.0d0
      do 10 k = 1,nd
         i = listd(k)
         ii = iky(i)+i
         e = e + 2.0d0*(h(ii)+f(ii))
         do 20 m = 1,ns
            j = lists(m)
            ij = itrian(i,j)
            e = e + g(int12(ij)+int34(ij))
 20      continue
 10   continue
c
      do 30 k = 1,ns
         i = lists(k)
         ii = iky(i) + i
         e = e + h(ii) + f(ii)
 30   continue
c
c     put this into the diagonal of hii
c
      e = 0.5d0*e + repnuc
      call dfill(numf,e,hii,numf1)
c
c     now add on the spin function dependent terms
c
      if (ns.ge.2 .and. nf(ns-2).gt.0) then
         do 40 k = 2,ns
            i = lists(k)
            do 50 m = 1,k-1
               j = lists(m)
               ij = itrian(i,j)
               gijij = g(int12(ij)+int34(ij))
               call selci_eijeji(work1,odonly,i,j,ns,indbar,iocc,
     $              w1,w2,work2)
               if (odonly) then
                  do kkk = 1, numfsq, numf1
                     hii(kkk) = hii(kkk) + gijij*work1(kkk)
                  enddo
               else
                  do kkk = 1, numfsq
                     hii(kkk) = hii(kkk) + gijij*work1(kkk)
                  enddo
               endif
 50         continue
 40      continue
      endif
c
      end
C>
C> @}
