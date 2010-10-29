      program main
*
* $Id$
*
      implicit none
c
c     The definitive MD5 test program is in mddriver.c
c     ... this is only to check the fortran interface 
c
      character*255 cdata
      integer n
      parameter (n=3350)
      double precision ddata(n)
      character*32 sum, sum1, sum10
      integer i
c
      cdata = 'This is so fulfilling ... '
      do i = 1, n
	 ddata(i) = i
      enddo
c
c     Compare simple with 1/10-at-a-time 
c
      call checksum_simple(8*n, ddata, sum)
      call checksum_init
      do i = 1, n
         call checksum_update(8, ddata(i))
      enddo
      call checksum_final(sum1)
      call checksum_init
      do i = 1, n, 10
         call checksum_update(80, ddata(i))
      enddo
      call checksum_final(sum10)
      write(6,*) ' simple ', sum
      write(6,*) ' 1x     ', sum1
      write(6,*) ' 10x    ', sum10
c
c     Check that changing 1 bit upsets things
c
      ddata(n) = ddata(n) + 1
      call checksum_simple(8*n, ddata, sum)
      write(6,*) ' diff1  ', sum
c
c     Check that changing order upsets things
c
      ddata(1) = 2
      ddata(2) = 1
      call checksum_simple(8*n, ddata, sum)
      write(6,*) ' order  ', sum
c
c     Check character interface
c
      call checksum_char_simple(cdata, sum)
      call checksum_init
      do i = 1, 255
         call checksum_char_update(cdata(i:i))
      enddo
      call checksum_final(sum1)
      call checksum_init
      do i = 1, 255, 5
         call checksum_char_update(cdata(i:i+4))
      enddo
      call checksum_final(sum10)
      write(6,*) ' simple ', sum
      write(6,*) ' 1x     ', sum1
      write(6,*) ' 5x     ', sum10
c
      end
