      double precision function util_random(iseed)
      implicit none
      integer iseed
#include "errquit.fh"
c $Id$
c
c routine to give a random number [0,1]
c
c     iseed is an integer seed when non-zero used to initialize the sequence     
c
#if defined(IBM) || defined(SP1) || ( defined(BGQ) && !defined(GFORTRAN) )
      if(iseed.gt.0) call srand(dble(iseed))
      util_random=rand()
#elif defined(CRAY)
      if(iseed.gt.0) call ranset(iseed)
      util_random=ranf()
#elif defined(LINUX)  || defined (MACX)
      double precision linux_rand
      external linux_rand
      if(iseed.gt.0) call linux_sran(iseed)
      util_random=linux_rand()
#else
      call errquit(' random routine missing', 0,CAPMIS_ERR)
#endif
      end


