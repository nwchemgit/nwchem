      double precision function util_random(iseed)
      implicit none
#include "errquit.fh"
c $Id: util_random.F,v 1.21 2003-10-17 22:58:35 carlfahl Exp $
c
c routine to give a random number [0,1]
c
c     iseed is an integer seed when non-zero used to initialize the sequence     
c
#if defined(IBM) || defined(SP1)
#define SET_RANDOM
      integer iseed
      if(iseed.gt.0) call srand(dble(iseed))
      util_random=rand()
#elif defined(KSR)
#define SET_RANDOM
      real random
      external random
      integer iseed
* on KSR if iseed != 0 it resets the sequence
      util_random=dble/work/u13377/nwchem(random(iseed))
#elif defined(SGI) || defined(SGITFP)
#define SET_RANDOM
       integer iseed
       real*8 rand
       external rand
       if(iseed.gt.0) call srand(dble(iseed))
       util_random=rand()
#elif defined(CRAY)
#define SET_RANDOM
      integer iseed
      if(iseed.gt.0) call ranset(iseed)
      util_random=ranf()
#elif defined(WIN32)
#define SET_RANDOM
      integer iseed
      double precision drand
      external drand
      if (iseed.gt.0) call srand(iseed)
      util_random = drand(0)
#elif (defined(LINUX) || defined(CYGNUS)) && !defined(LINUXALPHA)
#define SET_RANDOM
      integer iseed
      double precision linux_rand
      external linux_rand
      if(iseed.gt.0) call linux_sran(iseed)
      util_random=linux_rand()
#else
#define SET_RANDOM
      integer iseed
      double precision drand
      external drand
      util_random = drand(iseed)
#endif
#if !defined(SET_RANDOM)
      integer iseed
      call errquit('util_random not defined',0, UNKNOWN_ERR)
#endif
c
      end


