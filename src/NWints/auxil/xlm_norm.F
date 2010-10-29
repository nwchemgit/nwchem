      subroutine nwints_xlm_norm(qn,e,D,Ld)
* $Id$
      implicit none
#include "mafdecls.fh"
c
c     initialize data to normalize solid spherical harmonics
c
c     e(m) = -1 if m > 0 and m odd
c             1 otherwise
c     qn(m,l) = normalization constant required by Condon & Shortley
c               so that the unnormalized xml satisfy
c               qn( m,l) * x( m,l) = r^l * Re (Yml) m >= 0
c               qn(-m,l) * x(-m,l) = r^l * Im (Yml) m >= 0
c               |m| <= l, 0 otherwise
c
      integer Ld                       ! max angular momentum
      double precision e(-Ld:Ld)       ! [scratch] prefactor
      double precision qn(-Ld:Ld,0:Ld) ! [scratch] normalization constants
c......................................! [input/output] transformation matrix
      double precision D((((Ld+1)*(Ld+2))/2),-Ld:Ld,0:Ld) 
c::local
      integer L2s
      integer l, m
      double precision pi4
      double precision pi, fact
#if defined (QN_NEW_WAY)
      double precision pi, pi2x, resultlmm, resultlpm
#endif
      double precision norm

c
c... initialize e and qn

*rak:      write(6,*)' ma 1'
*rak:      if (.not.ma_verify_allocator_stuff()) stop ' ma ver'
      call dcopy((2*Ld+1)*(Ld+1),0.0d00,0,qn,1)
*rak:      write(6,*)' ma 2'
*rak:      if (.not.ma_verify_allocator_stuff()) stop ' ma ver'
*rak:      write(6,*) Ld
      do m = -Ld, 0
         e(m) = 1.0d0
      end do
      do m = 1, Ld
        if (mod(m,2).eq.0) then
          e(m) = 1.0d00
        else
          e(m) = -1.0d00
        endif
      end do

*rak:      write(6,*)' ma 3'
*rak:      if (.not.ma_verify_allocator_stuff()) stop ' ma ver'
c... form qn values via recursion

      pi4 = 16.0d0*atan(1.0d0)
      do l = 0, Ld
         qn(0,l) = (2*l+1)/pi4
         do m = 1, l
            qn( m,l) = qn(m-1,l)/((l+m)*(l-m+1))
            qn(-m,l) = qn(m,l)
         end do
      end do
*rak:      write(6,*)' ma 4'
*rak:      if (.not.ma_verify_allocator_stuff()) stop ' ma ver'
      do l = 0, Ld
         do m = -l, l
            qn(m,l) = e(m) * sqrt(qn(m,l))
         end do
      end do
*...  fix m=0 dependence (found emperically)
*...  multiply all m=0 coeffs by 1/sqrt(2)
*...  this gives all 2l+1 values the same diagonal value
      do l = 0, Ld
        qn(0,l) = qn(0,l)/sqrt(2.0d00)
      enddo
*... fix L>2 dependence (found emperically)
*...  multiply all coeffs by sqrt(8pi/(2l+1))
*...  this gives a diagional spherical overlap matrix
      pi = 2.0d0*dacos(0.0d00)
      do l = 0, Ld
        fact = 8.0d00*pi/(2.0d00*dble(l)+1.00d00)
        fact = sqrt(fact)    ! 
        do m = -l,l
          qn(m,l)=qn(m,l)*fact
        enddo
      enddo
#if defined (QN_NEW_WAY)
      call print_qn(Ld,qn,' old way ')
      pi = 2.0d0*dacos(0.0d00)
      pi2x = 2*pi
      write(6,*) 'pi =',pi
      call dcopy((2*Ld+1)*(Ld+1),0.0d00,0,qn,1)
      do l = 0, Ld
        do m = -l,l
          call dofact(resultlmm,(l-abs(m)))
          call dofact(resultlpm,(l+abs(m)))
          qn(m,l) = (2*l+1)/2 * resultlmm/resultlpm/pi2x
          qn(m,l) = e(m)*sqrt(qn(m,l))
        enddo
      enddo
      do l = 0, Ld
        qn(0,l) = qn(0,l)/sqrt(2.0d00)
      enddo
      call print_qn(Ld,qn,' new way ')
#endif
*rak:      write(6,*)' ma 5'
*rak:      if (.not.ma_verify_allocator_stuff()) stop ' ma ver'

c... transform D to normalized spherical coefficients

      do l = 0,Ld
        L2s = ((L+1)*(L+2))/2
        do m = -l,l
          norm = qn(m,l)
          call dscal(L2s,norm,D(1,m,l),1)
        enddo
      enddo
*rak:      write(6,*)' ma 6'
*rak:      if (.not.ma_verify_allocator_stuff()) stop ' ma ver'
c
      end
#if defined (QN_NEW_WAY)
      subroutine print_qn(Ld,qn,msg)
      implicit none
      integer Ld
      double precision qn(-Ld:Ld,0:Ld) 
      character*(*) msg

c
      integer l,m
      do l = 0, Ld
        do m = -l,l
          write(6,10000)msg,l,m,qn(m,l)

        enddo
      enddo
10000 format(1x,a,1x,'qn(',i3,',',i4,') = ',f20.6)
      end
      subroutine dofact(result,n)
      implicit none
      integer n
      double precision result
c
      integer i, nr
      nr = 1
      do i=n,1,-1
        nr = nr*i
      enddo
      result = dble(nr)
      write(6,*)' factorial of ',n,' is ', nr
      end
#endif
