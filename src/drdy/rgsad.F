C $Id: rgsad.F,v 1.2 2001-01-08 15:52:00 windus Exp $
c
c***********************************************************************
c     rgsad
c***********************************************************************
c
      subroutine drdy_rgsad(string,istrt,jtype)
c
c     Generic subroutine to read in the initial geometry for the saddle or 
c     starting stationary point. 
c
c     Jtype = 5  saddle point or starting stationary point
c
      implicit none             ! double precision (a-h,o-z)
#include "stdio.fh"
#include "drdyP.fh"
*::functions:
      integer  drdy_icint
      external drdy_icint
      double precision drdy_cfloat
      external         drdy_cfloat
*::passed:
      character * 80 string
      integer istrt
      integer jtype
*::local:
      integer i, isect, iend, ierr, ibegin
      integer j, jatm
      integer fu6
      fu6 = luout
c
c     read in one line of the geometry until natom have been read in
c
c
      do 1 i=1,natom
c
        call drdy_rline(string,istrt,isect,iend)
c
c        check for input error
c
        if (isect.eq.1.or.iend.eq.1) then
          write(fu6,1000)
1000      format(3x,'ERROR:  saddle point geometry specification',
     *          ' must end with and END')
          call errquit('drdy_rgsad: fatal error',911)
        end if
c
c        parse line to get atom number
c
        ibegin = istrt
        call drdy_rword(string,istrt,ierr)
        iend = istrt - 1
        jatm = drdy_icint(string(ibegin:iend))
c
c        check that jatm is a legal atom number
c
        if (jatm.le.0.or.jatm.gt.natom) then
          write(fu6,1100) string
1100      format(3x,'ERROR:  invalid atom number in saddle',
     *          ' point geometry specification.  Invalid line is: ',
     *          /1x,a80)
          call errquit('drdy_rgsad: fatal error',911)
        end if
c 
c        parse the x, y and z coordinates off the rest of the line
c
        do 2 j=2,0,-1
          ibegin = istrt
          call drdy_rword(string,istrt,ierr)
          iend = istrt - 1
c
c             check that the coordinate was found on the line
c
          if (ierr.eq.1.and.j.ne.0) then
            write(fu6,1200) string
1200        format(3x,'ERROR:  invalid cartesian coordinate in',
     *            ' saddle point.  Invalid line is: ',/a80)
            call errquit('drdy_rgsad: fatal error',911)
          end if
c
c            save coordinate
c
          xr(jatm*3-j,5) = drdy_cfloat(string(ibegin:iend))
c
2       continue
c
1     continue
      call drdy_rline(string,istrt,isect,iend)
      if(string(istrt:istrt+2).ne.'end') then
        write(fu6,1300)
        call errquit('drdy_rgsad: fatal error',911)
      end if
c
1300  format(2x,'ERROR: all atoms in the system must be given in the',
     *      ' starting point geometry specification')   
c
      return
      end
