      double precision function util_random(iseed)
      implicit none
c $Id: util_random.F,v 1.4 1997-04-18 23:03:10 d3e129 Exp $
c
c routine to give a random number [0,1]
c
c     iseed is an integer seed when non-zero used to initialize the sequence     
c
#if defined(IBM) || defined(SP1)
      integer iseed
      if(iseed.gt.0) call srand(dble(iseed))
      util_random=rand()
#elif defined(KSR)
      real random
      external random
      integer iseed
* on KSR if iseed != 0 it resets the sequence
      util_random=dble(random(iseed))
#elif defined(SGI) || defined(SGITFP)
      integer iseed
      real*8 rand
      external rand
      if(iseed.gt.0) call srand(dble(iseed))
      util_random=rand()
#elif defined(SUN) || defined(SOLARIS)
      integer iseed
      real*8 rand
      external rand
c      if(iseed.gt.0) call srand(dble(iseed))
      util_random=rand()
c #elif defined(__PARAGON__) || defined(DELTA)
c       integer iseed
c       if(iseed.gt.0) call srand(dble(iseed))
c       util_random=rand()
#elif defined(CRAY_T3D) || defined(CRAY_T3E)
      integer iseed
      if(iseed.gt.0) call ranset(iseed)
      util_random=ranf()
#else
      call errquit('util_random not defined',0)
#endif
c
      end


