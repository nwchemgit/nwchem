c
c $Id: qmmm_utils.F,v 1.4 2005-02-10 05:53:32 marat Exp $
c

      function qmmm_geom_distance(r1,r2)

      implicit none
      double precision r1(3),r2(3)
      double precision qmmm_geom_distance
c 
      integer k
      double precision r

      r=0.0d0
      do k=1,3
        r=r+(r1(k)-r2(k))*(r1(k)-r2(k))
      end do
      r=sqrt(r)

      qmmm_geom_distance=r

      end

      subroutine qmmm_sort(n,a)
      implicit none
      integer n
      integer a(n)
c
c     local variables:
      integer i
      integer pass  
      integer sorted 
      integer temp
      character*32 pname
      
      pass = 1
      sorted = 0
      do while(sorted .eq. 0) 
        sorted = 1
        do 2 i = 1,n-pass
          if(a(i) .gt. a(i+1)) then
            temp = a(i)
            a(i) = a(i+1)
            a(i+1) = temp
            sorted = 0
          endif
 2      continue
        pass = pass +1
      end do
      do i=1,n-1
       if(a(i).eq.a(i+1)) a(i)=-1 
      end do

      return
 
      end

      function qmmm_map(nr,ar,
     >                     n,a,map)
      implicit none
      integer nr
      integer ar(nr)
      integer n
      integer a(n)
      integer map(n)
      logical qmmm_map
c
c     local variables:
      integer i
      integer ir
      integer im

      qmmm_map = .true.
      im = 0
      do ir=1,nr
        do i=1,n
          if(a(i).eq.ar(ir)) then
            map(i)=ir
            im = im + 1 
          end if
        end do 
        if(im.eq.n) return
      end do

      qmmm_map = .false.
      return
      end

