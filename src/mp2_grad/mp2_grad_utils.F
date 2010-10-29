      subroutine mp2_zero_forces(force,oactive,nat)
*
* $Id$
*
      implicit none
      integer nat,i
      double precision force(3,nat)
      logical oactive
      dimension oactive(nat)
      do i=1,nat
         if (.not.oactive(i)) then
            force(1,i)=0.0d0
            force(2,i)=0.0d0
            force(3,i)=0.0d0
         endif
      enddo
      return
      end
c      
c$$$      subroutine div_eigen_aaaa(noa,nva,nob,
c$$$     $     nvb,nbf,orbe_a,t1_aaaa)
c$$$      implicit none
c$$$      integer noa,nva,nob,nvb,nbf
c$$$      double precision orbe_a(nbf)
c$$$      double precision t1_aaaa(*)
c$$$      integer a,a_count,b,b_count,count
c$$$c...  
c$$$c...  routine to return amplitude of the
c$$$c...  first order correction to the wavefunction for the
c$$$c...  alpha/alpha and beta/beta case.
c$$$c...  
c$$$      count=0
c$$$      a_count=0
c$$$      b_count=0
c$$$      do a=noa+1,nbf
c$$$         a_count=a_count+1
c$$$         do b=nob+1,nbf
c$$$            b_count=b_count+1
c$$$            count=count+1
c$$$            t1_aaaa(count)=t1_aaaa(count)/
c$$$     $           (orbe_a(i)+orbe_a(j)
c$$$     $           -orbe_a(a_count)-orbe_a(b_count))
c$$$         enddo
c$$$      enddo
c$$$      return
c$$$      end
c     
c     
      subroutine div_eigen_aabb(i,a_count,jlo_1,jhi_1,blo_1,bhi_1,
     $     noa,nva,nob,nvb,nbf,orbe_a,orbe_b,g_t1_aabb,t1_aabb)
      implicit none
#include "global.fh"
      integer noa,nva,nob,nvb,nbf
      double precision orbe_a(nbf),orbe_b(nbf)
      double precision t1_aabb(*)
      integer i,a_count
      integer j,j_count,b,b_count,count
      integer g_t1_aabb
      integer myid,jlo,jhi,blo,bhi
      integer blo_1,bhi_1,jlo_1,jhi_1
c...  
c...  routine to return amplitude of the
c...  first order correction to the wavefunction for the
c...  alpha/beta and beta/alpha case.
c...  
      myid=ga_nodeid()
      write(6,*) 'myid is ', myid
      write(6,*) 'INSIDE THE DIV_EIGEN_AABB CRAP'
      write(6,*) 'i,a_count,jlo_1,jhi_1,blo_1,bhi_1',
     $     i,a_count,jlo_1,jhi_1,blo_1,bhi_1
      call ga_print(g_t1_aabb)
c      call ga_distribution(g_t1_bb,myid,jlo,jhi,blo,bhi)
      count=0
      do j=jlo_1,jhi_1
         do b=blo_1,bhi_1
            count=count+1
            t1_aabb(count)=t1_aabb(count)/
     $           (orbe_a(i)+orbe_b(j)
     $           -orbe_a(a_count)-orbe_b(b))
         enddo
      enddo
      return
      end

