      subroutine tce_sort_0(unsorted,sorted,factor)
c
c $Id: tce_sort0.F,v 1.1 2008-09-25 14:54:18 jhammond Exp $
c
      implicit none
      double precision sorted
      double precision unsorted
      double precision factor
      sorted = unsorted * factor
      return
      end
