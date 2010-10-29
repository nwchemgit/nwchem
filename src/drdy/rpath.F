C $Id$
      subroutine drdy_rpath(string,iend,istrt)
c     Subroutine to read in the input for following the reaction path
      Implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "drdyP.fh"
*::functions:
      Double Precision drdy_cfloat
      External drdy_cfloat
*::passed
      Integer iend,istrt
      character * 80 string
*::local
      Integer ierr,isect
      character * 7 psign
      character * 5 minte
c     read in first keyword in this section
      call drdy_rline(string,istrt,isect,iend)
      do while (isect.eq.0.and.iend.eq.0)
c SCALEMASS
         if (string(istrt:istrt+8).eq.'scalemass') then
            call drdy_rword(string,istrt,ierr)
            if(ierr.ne.1) then
               redm = drdy_cfloat(string(istrt:80))
            else
               write(fu6,*)' ERROR:  variable SCALEMASS must ',
     *                    'have an argument'
               call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            end if
c SSTEP
         else if (string(istrt:istrt+4).eq.'sstep') then
            call drdy_rword(string,istrt,ierr)
            if (ierr.ne.1) then
               dels = drdy_cfloat(string(istrt:80))
            else
               write(fu6,*)' ERROR:  variable SSTEP must ',
     *                    'have an argument'
               call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            end if
c SSAVE
         else if (string(istrt:istrt+4).eq.'ssave') then
            call drdy_rword(string,istrt,ierr)
            if(ierr.ne.1) then
               delsv = drdy_cfloat(string(istrt:80))
            else
               write(fu6,*)' ERROR:  variable SSAVE must ',
     *                    'have an argument'
               call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            end if
c SHESS
         else if (string(istrt:istrt+4).eq.'shess') then
            call drdy_rword(string,istrt,ierr)
            if(ierr.ne.1) then
               delhss = drdy_cfloat(string(istrt:80))
            else
               write(fu6,*)' ERROR:  variable SHESS must ',
     *                    'have an argument'
               call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            end if
c SLP 
         else if (string(istrt:istrt+2).eq.'slp') then
            call drdy_rword(string,istrt,ierr)
            if(ierr.ne.1) then
               slp = drdy_cfloat(string(istrt:80))
            else
               write(fu6,*)' ERROR:  variable SLP must ',
     *                    'have an argument'
               call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            endif
c SLM 
         else if (string(istrt:istrt+2).eq.'slm') then
            call drdy_rword(string,istrt,ierr)
            if(ierr.ne.1) then
               slm = drdy_cfloat(string(istrt:80))
            else
               write(fu6,*)' ERROR:  variable SLM must ',
     *                    'have an argument'
               call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            endif
c SADDLE         
         else if (string(istrt:istrt+5).eq.'saddle') then   
            lgs(1) = 1
c NOSADDLE           
         else if (string(istrt:istrt+7).eq.'nosaddle') then    
            lgs(1) = 0
c SIGN
         else if (string(istrt:istrt+3).eq.'sign') then
            call drdy_rword(string,istrt,ierr)
            if(ierr.ne.1) then
               psign = string(istrt:istrt+6) 
               if (psign.eq.'reactan') lgs(27) = 0
               if (psign.eq.'product') lgs(27) = 1
            else
               write(fu6,*)' ERROR:  variable SIGN must ',
     *                    'have a character argument'
               call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            endif
c INTEGRA 
         else if (string(istrt:istrt+6).eq.'integra') then
            call drdy_rword(string,istrt,ierr)
            if (ierr.ne.1) then
              minte = string(istrt:istrt+4)
              if (minte.eq.'euler') lgs(31) = 1
              if (minte.eq.'lqa') lgs(31) = 2
              if (minte.eq.'clqa') lgs(31) = 3
              if (minte.eq.'cube') lgs(31) = 4
            else
              write (fu6,*) ' ERROR: variable INTEGRA must ',
     *              'have a character argument'
              call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
            endif
c PRINTFREQ
         else if (string(istrt:istrt+8).eq.'printfreq') then
            lgs(2) = 1
         else
            write(fu6,1000) string(istrt:80)
            call errquit('drdy_rpath: fatal error',911, INPUT_ERR)
         end if
         call drdy_rline(string,istrt,isect,iend)
      end do
c
1000  format(3x,'Error:  the following is not a valid keyword in the',
     *       ' path section',/A80)
c
      return
      end
