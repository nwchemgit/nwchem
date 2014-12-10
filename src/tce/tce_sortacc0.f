      subroutine tce_sortacc_0(unsorted,sorted,factor)
      implicit none
      double precision sorted
      double precision unsorted
      double precision factor
      sorted = sorted + unsorted * factor
      return
      end
