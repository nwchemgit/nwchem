c
c$Id$
      subroutine drdy_rstat(string,iend,istrt,jtype)
c     Generic subroutine to read in the input for all the stationary points.
c     Jtype = 1  first reactant 
c     Jtype = 2  second reactant 
c     Jtype = 3  first product 
c     Jtype = 4  second product
c     Jtype = 5  saddle point or starting stationary point
c
      Implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "drdyP.fh"
      Integer iend,istrt,jtype
      character * 80 string
C
      Integer ierr,isect
      character * 8  idmn 
c
c     read in first keyword in this section
c
      call drdy_rline(string,istrt,isect,iend)
c
      do while (isect.eq.0.and.iend.eq.0)
c GEOM
         if (string(istrt:istrt+3).eq.'geom') then
            if (jtype.le.4) call drdy_rgeom(string,istrt,jtype)
            if (jtype.eq.5) call drdy_rgsad(string,istrt,jtype)
c SPECIES
         else if (string(istrt:istrt+6).eq.'species') then              !0603PF97
            call drdy_rword(string,istrt,ierr)
            if (ierr.ne.1) then 
               idmn = string(istrt:istrt+7)
               if (idmn.eq.'atomic  ') icode(jtype)=1
               if (idmn.eq.'lints   '.or.idmn.eq.'linrp   ')
     >                icode(jtype)=3
               if (idmn.eq.'nonlints'.or.idmn.eq.'nonlinrp')
     >                icode(jtype) = 4
            else
               write(fu6,*) 'ERROR: variable SPECIES must have an',
     *                       ' argument!'
               call errquit('drdy_rstat: fatal error ',911, INPUT_ERR)
            end if
         else
            write(fu6,1000) string(istrt:80)
            call errquit('drdy_rstat: fatal error ',911, INPUT_ERR)
         endif
         call drdy_rline(string,istrt,isect,iend)
      end do
c
1000  format(3x,'Error:  the following is not a valid keyword in the',
     *       ' stationary points section',/A80)
c
      return
      end
