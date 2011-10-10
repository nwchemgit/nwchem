
      logical function is_sORd(l,lmax,locp)
      implicit none
      integer l
      integer lmax
      integer locp

      logical value
      integer d_start,d_end

      value = .false.

*     **** is an s function ****
      if ((locp.ne.0).and.(l.eq.1)) value = .true.

*     **** is an d function ****
      if ((locp.ne.2).and.(lmax.ge.2)) then

*        **** get the space spanned by d ****
         d_start = 5
         if (locp.eq.0) d_start=d_start-1
         if (locp.eq.1) d_start=d_start-3
         d_end = d_start + 5

*        **** check l against space spanned by d ****
         if ((l.ge.d_start).and.(l.le.d_end)) value = .true.
      end if

      is_sORd = value
      return
      end


c $Id$
