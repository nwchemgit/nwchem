      subroutine errquit(string, icode)
C$Id: errquit.F,v 1.6 1995-02-02 18:09:22 d3g681 Exp $
      character*(*) string
      integer icode
c     
c     error termination 
c     
c     string = error message printed to stdout
c     icode  = informative value printed to stdout
c     
      write(6,1)
 1    format(1x,72('-'))
      write(6, 2) string, icode
 2    format(1x,a,i4)
      write(6,1)
      call inp_errout
      write(6,1)
      write(6,*) ' current input line : '
      call inp_outrec
      write(6,1)
      call ga_error(string, icode)
*
#ifdef TCGMSG
      call parerr(icode)
#endif
*
#if defined(PARAGON)
	Call Kill( 0, 9)
#elif defined(DELTA)
        call killproc(-1,0)
#elif defined(IPSC)
        call killcube(-1,-1)
#endif
*
#if defined(SUN) || defined(SGI)
      call abort
#else
      stop 1
#endif
*
      end
