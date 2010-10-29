C $Id$
      subroutine drdy_rgeom(string,istrt,jtype)
c     jtype = 1  first reactant 
c     jtype = 2  second reactant 
c     jtype = 3  first product 
c     jtype = 4  second product
c
      Implicit none
#include "errquit.fh"
#include "drdyP.fh"
#include "stdio.fh"
*::functions:
      Double Precision drdy_cfloat
      External drdy_cfloat
      Integer drdy_icint
      External drdy_icint
*::passed:
      Integer istrt,jtype
      character * 80 string
*::local:
      Integer ibegin,iend,ierr,isect,jatm,nratm,j
c
c     read in each line of the geometry until an 'end' is found
c
      call drdy_rline(string,istrt,isect,iend)
      nratm = 1
c
      do while(string(istrt:istrt+2).ne.'end'.and.nratm.le.natom)
c
c        check for input error
c
         if (isect.eq.1.or.iend.eq.1) then
            write(fu6,1000)
1000        format(3x,'ERROR:  stationary point geometry specification',
     *               ' must end with and END')
            call errquit('drdy_rgeom: fatal error',911, INPUT_ERR)
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
1100         format(3x,'ERROR:  invalid atom number in stationary',
     *       ' point geometry specification.  Invalid line is: ',
     *       /1x,a80)
             call errquit('drdy_rgeom: fatal error',911, INPUT_ERR)
         end if
c
c        save atom number
c
         iatsv(nratm,jtype) = jatm
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
1200            format(3x,'ERROR:  invalid cartesian coordinate in',
     *          ' stationary point.  Invalid line is: ',/a80)
                call errquit('drdy_rgeom: fatal error',911, GEOM_ERR)
             end if
c
c            save coordinate
c
             if (jtype.eq.1.or.jtype.eq.2) 
     *           xr(jatm*3-j,1) = drdy_cfloat(string(ibegin:iend))
             if (jtype.eq.3.or.jtype.eq.4) 
     *           xr(jatm*3-j,3) = drdy_cfloat(string(ibegin:iend))
c
2        continue
c
9000     nratm = nratm+1
         call drdy_rline(string,istrt,isect,iend)
      end do
      nratom(jtype) = nratm - 1
      return
      end
