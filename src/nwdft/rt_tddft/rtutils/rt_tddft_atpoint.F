C
C     rt_tddft_atpoint.F
C
C     Determines if we are at a check/output point.
C
C     If npoints is at least one, we always check/print at the first
C     time step, then evenly distrubuted after that.  This ensures we
C     plot the ground state density matrix, etc.
C
      logical function rt_tddft_atpoint (it, nt, npoints)

      implicit none

      
C     == Inputs ==
      integer, intent(in) :: it, nt, npoints

      
C     == Variables ==
      logical needcheck

      
      if (npoints.lt.1) then
         needcheck = .false.
      elseif (npoints.eq.1) then
         needcheck = (it.eq.1)
      else
         needcheck = (mod(it, max(nt/(npoints-1), 1)) .eq. 0)
     $        .or. (it.eq.1)
      endif

      rt_tddft_atpoint = needcheck
      
      end function

c $Id$
