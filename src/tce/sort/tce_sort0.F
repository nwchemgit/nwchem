      subroutine tce_sort_0(unsorted,sorted,factor)
c
c $Id$
c
      implicit none
      double precision sorted
      double precision unsorted
      double precision factor
      sorted = unsorted * factor
      return
      end
