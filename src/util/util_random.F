      double precision function util_random(iseed)
      implicit none
c $Id: util_random.F,v 1.15 2000-03-13 19:43:40 edo Exp $
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
#endif
#if defined(KSR)
#define SET_RANDOM
      real random
      external random
      integer iseed
* on KSR if iseed != 0 it resets the sequence
      util_random=dble/work/u13377/nwchem(random(iseed))
#endif
#if defined(SGI) || defined(SGITFP)
#define SET_RANDOM
       integer iseed
       real*8 rand
       external rand
       if(iseed.gt.0) call srand(dble(iseed))
       util_random=rand()
#endif
#if defined(HPUX)
#define SET_RANDOM
       integer iseed
       real*8 drand
       external drand
       util_random=drand()
#endif
#if defined(SUN) || defined(SOLARIS)
#define SET_RANDOM
      integer iseed
c      integer rand
c      external rand
      real*8 drand
      external drand
c      if(iseed.gt.0) call srand(iseed)
c      util_random=3.051850948d-05*dble(rand())
      util_random=drand(iseed)
#endif
c #elif defined(__PARAGON__) || defined(DELTA)
c       integer iseed
c       if(iseed.gt.0) call srand(dble(iseed))
c       util_random=rand()
#if defined(CRAY_T3D) || defined(CRAY_T3E)
#define SET_RANDOM
      integer iseed
      if(iseed.gt.0) call ranset(iseed)
      util_random=ranf()
#endif
#if (defined(LINUX) || defined(CYGNUS)) && !(defined(__linux__) && defined(__alpha)) 
#define SET_RANDOM
      integer iseed
      double precision linux_rand
      external linux_rand
      if(iseed.gt.0) call linux_sran(iseed)
      util_random=linux_rand()
#endif
#if defined(DECOSF) || (defined(__linux__) && defined(__alpha)) 
#define SET_RANDOM
      integer iseed
      double precision drand
      external drand
      if (iseed.gt.0) call srand(iseed)
      util_random = drand()
#endif
#if defined(WIN32)
#define SET_RANDOM
      integer iseed
      double precision drand
      external drand
      if (iseed.gt.0) call srand(iseed)
      util_random = drand(0)
#endif
#if !defined(SET_RANDOM)
      integer iseed
      call errquit('util_random not defined',0)
#endif
c
      end


