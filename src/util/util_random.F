      double precision function util_random(iseed)
      implicit none
c $Id: util_random.F,v 1.1 1997-04-16 18:31:11 d3j191 Exp $
c
c routine to give a random number [0,1]
c
c     iseed is an integer seed when non-zero used to initialize the sequence     
c
#if defined(IBM) || defined(SP1)
      integer iseed
      if(iseed.gt.0) call srand(dble(iseed))
      util_random=rand()
c #elif defined(KSR)
c       integer*8 iseed
c       util_random=random(iseed)
#elif defined(SUN) || defined(SGI) || defined(SGITFP) || defined(SOLARIS)
      integer iseed
      real*8 rand
      external rand
      if(iseed.gt.0) call srand(dble(iseed))
      util_random=rand()
c #elif defined(__PARAGON__) || defined(DELTA)
c       integer iseed
c       if(iseed.gt.0) call srand(dble(iseed))
c       util_random=rand()
c #elif defined(CRAY_T3D) || defined(CRAY_T3E)
c       integer iseed
c       if(iseed.gt.0) call srand(dble(iseed))
c       util_random=rand()
#else
      call errquit('util_random not defined',0)
#endif
c
      end


