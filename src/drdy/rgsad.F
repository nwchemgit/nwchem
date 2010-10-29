C $Id$
      subroutine drdy_rgsad(string,istrt,jtype)
c
c     Generic subroutine to read in the initial geometry for the saddle or 
c     starting stationary point. 
c
c     Jtype = 5  saddle point or starting stationary point
c
      Implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "drdyP.fh"
*::functions:
      Integer  drdy_icint
      External drdy_icint
      Double Precision drdy_cfloat
      External         drdy_cfloat
*::passed:
      Integer istrt,jtype
      character * 80 string
*::local:
      Integer i,ibegin,iend,ierr,isect,j,jatm
c
c     read in one line of the geometry until natom have been read in
c
c
      do i=1,natom
         call drdy_rline(string,istrt,isect,iend)
c
c        check for input error
c
         if (isect.eq.1.or.iend.eq.1) then
            write(fu6,1000)
            call errquit('drdy_rgsad: fatal error',911, INPUT_ERR)
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
             call errquit('drdy_rgsad: fatal error',911, INPUT_ERR)
         end if
c 
c        parse the x, y and z coordinates off the rest of the line
c
         do j=2,0,-1
             ibegin = istrt
             call drdy_rword(string,istrt,ierr)
             iend = istrt - 1
c
c             check that the coordinate was found on the line
c
             if (ierr.eq.1.and.j.ne.0) then
                write(fu6,1200) string
                call errquit('drdy_rgsad: fatal error',911, INPUT_ERR)
             end if
c
c            save coordinate
c
             xr(jatm*3-j,5) = drdy_cfloat(string(ibegin:iend))
c
          enddo
c
      enddo
      call drdy_rline(string,istrt,isect,iend)
      if(string(istrt:istrt+2).ne.'end') then
        write(fu6,1300)
        call errquit('drdy_rgsad: fatal error',911, INPUT_ERR)
      end if
c
      return
1000        format(3x,'ERROR:  saddle point geometry specification',
     *               ' must end with and END')
1100         format(3x,'ERROR:  invalid atom number in saddle',
     *       ' point geometry specification.  Invalid line is: ',
     *       /1x,a80)
1200            format(3x,'ERROR:  invalid cartesian coordinate in',
     *          ' saddle point.  Invalid line is: ',/a80)
1300  format(2x,'ERROR: all atoms in the system must be given in the',
     * ' starting point geometry specification')   
      end
