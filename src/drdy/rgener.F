C $Id$
      Subroutine drdy_rgener(string,iend,istrt)
c
c     called by:
c           read5
c     calls:
c           drdy_rline
c
      Implicit none 
#include "errquit.fh"
#include "drdyP.fh"
c           rword
C
      Integer iend,istrt,isect
      Character * 80 string
*::local:
      Integer ierr,lvec,lhess,lspc
c
c     read in first keyword in this section
c
      call drdy_rline(string,istrt,isect,iend)
c
      do while (isect.eq.0.and.iend.eq.0)
c TITLE
         if (string(istrt:istrt+4).eq.'title') then
*            READ(5,500) ctitle
            READ(lu_drdy_rline,500) ctitle
c RESTART
* should never appear in nwchem input deck
         else if (string(istrt:istrt+5).eq.'restar') then
              lgs(8) = 1
c RE-START
         else if (string(istrt:istrt+6).eq.'re-star') then
              lgs(8) = 1
c ATOMS
         else if (string(istrt:istrt+4).eq.'atoms') then
            call drdy_ratoms(string,istrt)
c SINGLE POINT CALCULATION
         else if (string(istrt:istrt+10).eq.'singlepoint') then
            lgs(39) = 1
C SAVE ELECTRONIC STRUCTURE FILES
         else if (string(istrt:istrt+7).eq.'savefile') then
            lvec = 0
            lhess = 0
            lspc = 0
            call drdy_rword(string,istrt,ierr)
            do while (ierr.eq.0)
               if(string(istrt:istrt+3).eq.'vecs') then
                  lvec = 1
               elseif(string(istrt:istrt+3).eq.'hess') then
                  lhess = 1
               elseif(string(istrt:istrt+2).eq.'spc') then
                  if(lgs(39).ne.0) then
                     lspc = 1
                  else
                     write(fu6,1202)
                  endif
               else
                  write(fu6,1201) string(istrt:istrt+2)
                  call errquit('drdy_rgener: invalid argument',555,
     &       INPUT_ERR)
               endif
               call drdy_rword(string,istrt,ierr)
            enddo
            if (lhess.eq.0) then
               if (lvec.eq.0) then
                  write(fu6,1203)
                  call errquit('drdy_rgener: invalid argument',555,
     &       INPUT_ERR)
               else
                  lgs(38) = 1
               endif
               if(lspc.ne.0) write(fu6,1204)
            else
               if (lvec.eq.0) then
                  lgs(38) = 2
                  if(lspc.ne.0) write(fu6,1204)
               else if (lspc.eq.0) then
                  lgs(38) = 3
               else
                  lgs(38) = 4
               endif
            endif
         else
            write(fu6,1200) string(istrt:80)
            call errquit('drdy_rgener: fatal error ',911, UNKNOWN_ERR)
         end if
         call drdy_rline(string,istrt,isect,iend)
      end do
c
500   format(a80)
1200  format(3x,'Error:  the following is not a valid keyword in the',
     *       ' general section',/A80)
1201  format(3x,'Error:  the following is not a valid argument for the',
     *       ' SAVEFILE keyword',/,A80)
1202  Format(3x,'Cannot save vectors for single point calculation ',
     *       'because the calculation is not being performed')
1203  format(3x,'SAVEFILE is specified, but no valid arguments are',
     *       ' given')
1204  format(3x,'To save the vectors from single point calculation,',
     *       ' the regular vectors and hessian files must also be',
     *       ' saved')
c
      return
      end
