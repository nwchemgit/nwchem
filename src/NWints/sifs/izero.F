      subroutine izero(n,ia,inca)
c
c zero the elements of ia
c
      integer n, inca
      integer ia(inca,*)
c
      integer i
      do 10 i = 1, n
         ia(1,i) = 0
10    continue
c
      return
      end
