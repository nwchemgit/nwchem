      block data xlm_block_data
      implicit none
#include "xlm.fh"
c
c     xlm_max25(m,l) = absolute maximum value of blm * xlm in a cube
c     (-1,+1) found by brute force search.  The maximum value is
c     always realized on the boundary.
c
      data xlm_max25(1)/
c     l = 0
     $     1.0000d+00/
      data xlm_max25(2:4)/
c     l = 1
     $     5.0000d-01, 1.0000d+00, 5.0000d-01/
      data xlm_max25(5:9)/
c     l = 2
     $     2.5000d-01, 5.0000d-01, 5.0000d-01, 5.0000d-01, 1.2500d-01/
      data xlm_max25(10:16)/
c     l = 3
     $     4.1667d-02, 2.5000d-01, 1.8750d-01, 3.3333d-01, 1.8750d-01,
     $     1.2500d-01,
     $     4.1667d-02/
      data xlm_max25(17:25)/
c     l = 4
     $     4.0094d-03, 4.1667d-02, 8.3333d-02, 5.8926d-02, 1.4583d-01,
     $     5.8926d-02,
     $     5.2083d-02, 4.1667d-02, 1.0417d-02/
      data xlm_max25(26:36)/
c     l = 5
     $     1.0417d-03, 4.0094d-03, 1.5886d-02, 1.6038d-02, 3.1250d-02,
     $     2.1108d-02,
     $     3.1250d-02, 1.0417d-02, 1.5886d-02, 1.0417d-02, 1.0417d-03/
      data xlm_max25(37:49)/
c     l = 6
     $     1.7361d-04, 1.0417d-03, 1.7391d-03, 2.7141d-03, 7.8125d-03,
     $     6.2500d-03,
     $     8.9104d-03, 6.2500d-03, 1.3459d-03, 2.7141d-03, 4.1667d-03,
     $     1.0417d-03,
     $     6.2618d-05/
      data xlm_max25(50:64)/
c     l = 7
     $     1.2401d-05, 1.7361d-04, 4.3489d-04, 4.0764d-04, 7.1615d-04,
     $     2.2569d-03,
     $     9.7131d-04, 2.9762d-03, 9.7131d-04, 4.6600d-04, 7.1615d-04,
     $     6.9444d-04,
     $     4.3489d-04, 6.2618d-05, 1.2401d-05/
      data xlm_max25(65:81)/
c     l = 8
     $     5.4256d-07, 1.2401d-05, 7.4405d-05, 9.0788d-05, 4.9359d-05,
     $     2.6042d-04,
     $     2.0488d-04, 3.5962d-04, 3.1603d-04, 3.5962d-04, 1.0353d-04,
     $     2.6042d-04,
     $     1.1285d-04, 9.0788d-05, 2.7902d-05, 1.2401d-05, 1.5501d-06/
      data xlm_max25(82:100)/
c     l = 9
     $     8.6117d-08, 5.4256d-07, 5.4280d-06, 1.6534d-05, 8.6967d-06,
     $     7.9191d-06,
     $     3.6873d-05, 5.7034d-05, 6.6654d-05, 6.3734d-05, 6.6654d-05,
     $     1.3176d-05,
     $     3.6873d-05, 4.3403d-05, 8.6967d-06, 7.0573d-06, 5.4280d-06,
     $     1.5501d-06,
     $     8.6117d-08/
      data xlm_max25(101:121)/
c     l = 10
     $     8.6117d-09, 8.6117d-08, 2.4666d-07, 1.3029d-06, 1.4519d-06,
     $     2.7459d-06,
     $     2.2321d-06, 3.9691d-06, 1.5243d-05, 5.0225d-06, 1.7912d-05,
     $     5.0225d-06,
     $     1.5765d-06, 3.9691d-06, 8.2672d-06, 2.7459d-06, 1.0330d-06,
     $     1.3029d-06,
     $     6.8893d-07, 8.6117d-08, 2.9625d-09/
      data xlm_max25(122:144)/
c     l = 11
     $     3.9144d-10, 8.6117d-09, 3.8759d-08, 6.5900d-08, 1.6741d-07,
     $     3.0426d-07,
     $     6.8463d-07, 3.1335d-07, 1.1626d-06, 1.4640d-06, 1.3633d-06,
     $     1.9729d-06,
     $     1.3633d-06, 3.0222d-07, 1.1626d-06, 5.2268d-07, 6.8463d-07,
     $     7.4962d-08,
     $     1.6741d-07, 1.7223d-07, 3.8759d-08, 2.9625d-09, 3.9144d-10/
      data xlm_max25(145:169)/
c     l = 12
     $     1.1095d-11, 3.9144d-10, 3.9144d-09, 1.0073d-08, 1.0947d-08,
     $     1.5285d-08,
     $     8.8987d-08, 6.9496d-08, 2.4395d-08, 1.7510d-07, 1.6743d-07,
     $     2.5444d-07,
     $     1.8335d-07, 2.5444d-07, 3.4102d-08, 1.7510d-07, 1.3095d-07,
     $     6.9496d-08,
     $     1.9492d-08, 1.5285d-08, 2.3682d-08, 1.0073d-08, 1.3666d-09,
     $     3.9144d-10,
     $     3.2620d-11/
      data xlm_max25(170:196)/
c     l = 13
     $     1.2546d-12, 1.1095d-11, 1.7942d-10, 1.0438d-09, 1.5681d-09,
     $     1.0419d-09,
     $     4.6647d-09, 1.1482d-08, 6.6611d-09, 4.2914d-09, 1.0560d-08,
     $     3.9666d-08,
     $     1.5797d-08, 4.3610d-08, 1.5797d-08, 2.8175d-09, 1.0560d-08,
     $     2.7629d-08,
     $     6.6611d-09, 3.7185d-09, 4.6647d-09, 1.4025d-09, 1.5681d-09,
     $     3.7933d-10,
     $     1.7942d-10, 3.2620d-11, 1.2546d-12/
      data xlm_max25(197:225)/
c     l = 14
     $     8.9615d-14, 1.2546d-12, 5.1743d-12, 4.8980d-11, 1.7125d-10,
     $     1.2767d-10,
     $     8.5919d-11, 7.6890d-10, 6.2152d-10, 1.6147d-09, 6.2564d-10,
     $     2.1018d-09,
     $     4.2324d-09, 2.4087d-09, 5.1188d-09, 2.4087d-09, 4.7113d-10,
     $     2.1018d-09,
     $     2.1529d-09, 1.6147d-09, 3.6558d-10, 7.6890d-10, 4.0324d-10,
     $     1.2767d-10,
     $     6.8493d-11, 4.8980d-11, 1.5055d-11, 1.2546d-12, 3.0241d-14/
      data xlm_max25(226:256)/
c     l = 15
     $     2.9872d-15, 8.9615d-14, 5.8252d-13, 1.4763d-12, 8.5307d-12,
     $     1.4679d-11,
     $     1.5959d-11, 2.1628d-11, 5.6684d-11, 1.4434d-10, 1.8838d-10,
     $     5.0945d-11,
     $     3.4659d-10, 2.4343d-10, 4.4319d-10, 2.6267d-10, 4.4319d-10,
     $     4.6605d-11,
     $     3.4659d-10, 1.9604d-10, 1.8838d-10, 2.2105d-11, 5.6684d-11,
     $     8.3880d-11,
     $     1.5959d-11, 7.8739d-12, 8.5307d-12, 4.1820d-12, 5.8252d-13,
     $     3.0241d-14,
     $     2.9872d-15/
      data xlm_max25(257:289)/
c     l = 16
     $     6.2632d-17, 2.9872d-15, 4.1820d-14, 1.6437d-13, 2.8208d-13,
     $     9.0383d-13,
     $     1.1759d-12, 3.8445d-12, 2.8582d-12, 5.7786d-12, 2.3533d-11,
     $     9.5102d-12,
     $     4.0292d-12, 2.0605d-11, 4.8337d-11, 3.0675d-11, 5.1133d-11,
     $     3.0675d-11,
     $     3.1937d-12, 2.0605d-11, 3.9552d-11, 9.5102d-12, 4.6078d-12,
     $     5.7786d-12,
     $     7.4269d-12, 3.8445d-12, 4.7901d-13, 9.0383d-13, 7.5426d-13,
     $     1.6437d-13,
     $     1.4221d-14, 2.9872d-15, 1.8670d-16/
      data xlm_max25(290:324)/
c     l = 17
     $     5.4911d-18, 6.2632d-17, 1.4003d-15, 1.1949d-14, 3.0712d-14,
     $     3.6509d-14,
     $     4.5191d-14, 2.9872d-13, 4.3441d-13, 2.0103d-13, 1.1504d-12,
     $     1.3980d-12,
     $     1.6549d-12, 6.2532d-13, 2.0666d-12, 5.7712d-12, 2.3341d-12,
     $     6.5798d-12,
     $     2.3341d-12, 4.1179d-13, 2.0666d-12, 3.7056d-12, 1.6549d-12,
     $     5.1498d-13,
     $     1.1504d-12, 4.3174d-13, 4.3441d-13, 6.5488d-14, 4.5191d-14,
     $     8.5134d-14,
     $     3.0712d-14, 4.1415d-15, 1.4003d-15, 1.8670d-16, 5.4911d-18/
      data xlm_max25(325:361)/
c     l = 18
     $     3.0506d-19, 5.4911d-18, 2.9649d-17, 4.0462d-16, 2.2871d-15,
     $     3.7971d-15,
     $     2.9369d-15, 9.8933d-15, 4.1634d-14, 2.3330d-14, 1.4976d-14,
     $     1.0196d-13,
     $     1.2401d-13, 2.3076d-13, 5.0792d-14, 3.4427d-13, 2.2079d-13,
     $     4.0429d-13,
     $     2.7181d-13, 4.0429d-13, 3.8895d-14, 3.4427d-13, 1.6576d-13,
     $     2.3076d-13,
     $     3.1796d-14, 1.0196d-13, 8.9615d-14, 2.3330d-14, 1.2280d-14,
     $     9.8933d-15,
     $     4.0921d-15, 3.7971d-15, 8.2386d-16, 4.0462d-16, 8.7858d-17,
     $     5.4911d-18,
     $     1.0188d-19/
      data xlm_max25(362:400)/
c     l = 19
     $     8.0280d-21, 3.0506d-19, 2.5931d-18, 8.7721d-18, 7.9282d-17,
     $     2.9561d-16,
     $     2.7437d-16, 1.3919d-16, 1.7114d-15, 2.4686d-15, 2.9581d-15,
     $     2.5382d-15,
     $     4.3652d-15, 2.1470d-14, 1.2176d-14, 2.6276d-15, 2.3757d-14,
     $     3.3066d-14,
     $     3.2371d-14, 3.4813d-14, 3.2371d-14, 2.3408d-15, 2.3757d-14,
     $     2.8558d-14,
     $     1.2176d-14, 2.7834d-15, 4.3652d-15, 1.0870d-14, 2.9581d-15,
     $     1.1711d-15,
     $     1.7114d-15, 6.0867d-16, 2.7437d-16, 1.1529d-16, 7.9282d-17,
     $     2.5625d-17,
     $     2.5931d-18, 1.0188d-19, 8.0280d-21/
      data xlm_max25(401:441)/
c     l = 20
     $     1.3357d-22, 8.0280d-21, 1.4450d-19, 7.6277d-19, 1.7981d-18,
     $     1.0783d-17,
     $     2.2880d-17, 1.7696d-17, 3.0589d-17, 1.3820d-16, 1.8423d-16,
     $     4.6675d-16,
     $     2.1894d-16, 7.6255d-16, 1.6803d-15, 1.0010d-15, 3.5723d-16,
     $     1.2297d-15,
     $     4.3069d-15, 1.3708d-15, 4.7251d-15, 1.3708d-15, 2.1405d-16,
     $     1.2297d-15,
     $     3.1622d-15, 1.0010d-15, 3.5874d-16, 7.6255d-16, 4.7147d-16,
     $     4.6675d-16,
     $     5.8691d-17, 1.3820d-16, 1.2252d-16, 1.7696d-17, 1.0977d-17,
     $     1.0783d-17,
     $     5.1098d-18, 7.6277d-19, 4.8484d-20, 8.0280d-21, 4.0140d-22/
      data xlm_max25(442:484)/
c     l = 21
     $     9.5571d-24, 1.3357d-22, 3.8133d-21, 4.2816d-20, 1.5466d-19,
     $     2.6574d-19,
     $     9.6749d-19, 9.7500d-19, 3.8323d-18, 3.9781d-18, 5.4271d-18,
     $     3.3191d-17,
     $     3.1525d-17, 1.0318d-17, 8.7052d-17, 6.6373d-17, 1.4819d-16,
     $     2.9862d-17,
     $     1.9078d-16, 1.9825d-16, 2.0937d-16, 2.3685d-16, 2.0937d-16,
     $     2.0565d-17,
     $     1.9078d-16, 1.1467d-16, 1.4819d-16, 2.3986d-17, 8.7052d-17,
     $     4.7694d-17,
     $     3.1525d-17, 6.1368d-18, 5.4271d-18, 1.1885d-17, 3.8323d-18,
     $     6.0016d-19,
     $     9.6749d-19, 7.1690d-19, 1.5466d-19, 1.4526d-20, 3.8133d-21,
     $     4.0140d-22,
     $     9.5571d-24/
      data xlm_max25(485:529)/
c     l = 22
     $     4.3441d-25, 9.5571d-24, 6.3845d-23, 1.1374d-21, 8.7973d-21,
     $     2.2422d-20,
     $     2.8150d-20, 4.5643d-20, 2.0809d-19, 4.4906d-19, 2.7905d-19,
     $     9.7432d-19,
     $     2.9712d-18, 1.4144d-18, 1.0955d-18, 3.9283d-18, 1.0227d-17,
     $     9.3680d-18,
     $     1.5061d-18, 1.5554d-17, 1.4450d-17, 1.9390d-17, 1.5144d-17,
     $     1.9390d-17,
     $     1.1581d-18, 1.5554d-17, 1.2610d-17, 9.3680d-18, 1.0615d-18,
     $     3.9283d-18,
     $     7.0135d-18, 1.4144d-18, 8.0976d-19, 9.7432d-19, 4.4321d-19,
     $     4.4906d-19,
     $     4.7605d-20, 4.5643d-20, 6.7970d-20, 2.2422d-20, 3.0470d-21,
     $     1.1374d-21,
     $     1.9114d-22, 9.5571d-24, 1.4413d-25/
      data xlm_max25(530:576)/
c     l = 23
     $     9.4438d-27, 4.3441d-25, 4.5614d-24, 1.9324d-23, 2.3665d-22,
     $     1.3045d-21,
     $     2.2906d-21, 1.9702d-21, 5.2336d-21, 2.9871d-20, 2.6399d-20,
     $     9.9482d-21,
     $     1.1587d-19, 9.3278d-20, 2.2458d-19, 1.1431d-19, 3.0305d-19,
     $     9.9745d-19,
     $     3.9380d-19, 1.2462d-19, 4.7541d-19, 1.9182d-18, 5.9511d-19,
     $     2.0446d-18,
     $     5.9511d-19, 6.9785d-20, 4.7541d-19, 1.5465d-18, 3.9380d-19,
     $     1.4462d-19,
     $     3.0305d-19, 4.4452d-19, 2.2458d-19, 5.4170d-20, 1.1587d-19,
     $     5.2011d-20,
     $     2.6399d-20, 8.4540d-21, 5.2336d-21, 3.3259d-21, 2.2906d-21,
     $     4.6878d-22,
     $     2.3665d-22, 5.7343d-23, 4.5614d-24, 1.4413d-25, 9.4438d-27/
      data xlm_max25(577:625)/
c     l = 24
     $     1.3023d-28, 9.4438d-27, 2.0776d-25, 1.3757d-24, 4.1284d-24,
     $     3.5907d-23,
     $     1.3858d-22, 1.4792d-22, 6.7697d-23, 8.8483d-22, 2.0596d-21,
     $     1.3235d-21,
     $     1.6382d-21, 6.0169d-21, 1.2072d-20, 2.0276d-20, 6.2653d-21,
     $     3.9495d-20,
     $     3.2449d-20, 5.4901d-20, 1.1198d-20, 6.4038d-20, 1.0950d-19,
     $     6.9376d-20,
     $     1.2389d-19, 6.9376d-20, 7.1508d-21, 6.4038d-20, 7.2756d-20,
     $     5.4901d-20,
     $     1.0512d-20, 3.9495d-20, 1.7108d-20, 2.0276d-20, 1.9322d-21,
     $     6.0169d-21,
     $     7.6396d-21, 1.3235d-21, 7.9166d-22, 8.8483d-22, 2.5129d-22,
     $     1.4792d-22,
     $     5.3387d-23, 3.5907d-23, 1.2055d-23, 1.3757d-24, 6.9148d-26,
     $     9.4438d-27,
     $     3.9349d-28/
      data xlm_max25(626:676)/
c     l = 25
     $     7.8698d-30, 1.3023d-28, 4.5252d-27, 6.2958d-26, 2.9206d-25,
     $     6.5497d-25,
     $     3.9917d-24, 9.5571d-24, 5.7688d-24, 1.1315d-23, 7.3927d-23,
     $     6.4324d-23,
     $     2.1462d-22, 1.6135d-22, 3.0917d-22, 1.4692d-21, 7.6723d-22,
     $     2.8935d-22,
     $     2.1781d-21, 3.0885d-21, 4.2552d-21, 5.6434d-22, 6.1233d-21,
     $     4.2730d-21,
     $     7.1961d-21, 4.4601d-21, 7.1961d-21, 3.9650d-22, 6.1233d-21,
     $     3.7682d-21,
     $     4.2552d-21, 4.4785d-22, 2.1781d-21, 2.4024d-21, 7.6723d-22,
     $     2.7039d-22,
     $     3.0917d-22, 5.0830d-22, 2.1462d-22, 4.0779d-23, 7.3927d-23,
     $     4.6790d-23,
     $     5.7688d-24, 4.3368d-24, 3.9917d-24, 1.8607d-24, 2.9206d-25,
     $     2.1107d-26,
     $     4.5252d-27, 3.9349d-28, 7.8698d-30/
c
      end
      subroutine xlm_norm(lmax, q, lq)
      implicit none
#include "xlm.fh"
      integer lmax, lq
      double precision q(-lq:lq,0:lmax)
c     
c     Return in q the normalization constants qn ... see xlm_init
c
      integer l, m
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = qn(m,l)
         enddo
      enddo
c
      end
      subroutine xlm_thresh(lmax, q, lq, thresh, dist, luse)
      implicit none
#include "xlm.fh"
      integer lmax, lq, luse
      double precision q(-lq:lq,0:lmax), thresh, dist
c
c     Return in luse the angular momentum required to obtain
c     the required accuracy in the potential at a distance dist.
c
c     For use in the FMM dist should be (ws+0.5)*box.  
c
c     If all components are negligble luse=-1 is returned.
c
      integer l, m, loff
      double precision scale, test
c
      scale = 4.0d0/(2.0d0*dist)**(lmax+1)
      do l = lmax,0,-1
         loff = l*l + l + 1
         do m = -l, l
            test = abs(scale*q(m,l)*c(m,l)*xlm_max25(loff+m))
            if (test .gt. thresh) goto 10
         enddo
         scale = scale * dist * 2.0d0
      enddo
 10   luse = l
c
      end
      subroutine xlm(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      integer lmax, lq
      double precision x, y, z, q(-lq:lq, 0:lmax)
c
      integer l, m, lp1
      double precision r2, twolp1
c
c     compute the unnormalized real solid spherical harmonics
c     up to order lmax at the point x, y, z.
c     lq specifies the dimension of q.
c
      if (lmax .gt. lq) then
         write(6,*) ' invalid dimension for q ', lmax, lq
         stop 1
      end if
***      call dfill((2*lq+1)*(lmax+1), 0.0d0, q, 1)
      do l = 0,lmax             ! NOTE that |m| = l+1 is referenced
         do m = -lq,lq
            q(m,l) = 0.0d0
         end do
      end do
      q(0,0) = 1.0d0            ! MUST be after zero
      if (lmax.ge.1) then
         q(-1,1) = y
         q( 0,1) = z
         q( 1,1) = x
      end if
c
      r2 = x*x + y*y + z*z
      do l = 1, lmax-1
         lp1 = l + 1
         twolp1 = 2*l + 1
         q( lp1, lp1) = twolp1*(x*q( l,l) - y*q(-l,l))
         q(-lp1, lp1) = twolp1*(x*q(-l,l) + y*q( l,l))
         q(0,lp1) = (twolp1*z*q(0,l) - l*r2*q(0,l-1))*fastdiv(l+1)
c$$$         do m = 1,l
c$$$            q( m,lp1) = (twolp1*z*q( m,l) - (l+m)*r2*q( m,l-1))/(l-m+1)
c$$$            q(-m,lp1) = (twolp1*z*q(-m,l) - (l+m)*r2*q(-m,l-1))/(l-m+1)
c$$$         end do
         do m = 1,l
            q( m,lp1) = (twolp1*z*q( m,l) - (l+m)*r2*q( m,l-1)) *
     $           fastdiv(l-m+1)
            q(-m,lp1) = (twolp1*z*q(-m,l) - (l+m)*r2*q(-m,l-1)) *
     $           fastdiv(l-m+1)
         end do
      end do
c
      end
      double precision function factorial(n)
      implicit none
      integer n, i
      double precision f
c
      f = 1.0d0
      do i = 1, n
         f = f * dble(i)
      enddo
c
      factorial = f
c
      end
      double precision function double_factorial(n)
      implicit none
      integer n, i
      double precision f
c
      f = 1.0d0
      do i = n, 1, -2
         f = f * dble(i)
      enddo
c
      double_factorial = f
c
      end
      subroutine xlm_coeff_inv(lmax,d,dinv)
      implicit none
      integer lmax
      double precision d(((lmax+1)*(lmax+2))/2, -lmax:lmax, 0:lmax)
      double precision 
     $     dinv(((lmax+1)*(lmax+2)*(lmax+3))/6,-lmax:lmax,0:lmax)
c
c     NOTE differing dimensions for d and dinv
c
c     Compute the coefficients of the solid harmonics for expanding
c     the cartesian basis ... the inverse of xlm_coeff.  I.e., return
c     in dinv
c
c     x^i y^j z^k = sum(l=i+j+k,0,-2) sum(m=-l..l) 
c     .                     r^(l-i-j-k) dinv(ijk,m,l) Xlm
c
c     Also returned (as an extra, added bonus for the wary shopper)
c     is d from xlm_coeff.
c
c     Xlm = sum(ijk) x^i y^j z^k d(ijk,m,l)
c     
c     dinv(ijk,lm) = <Xlm * x^iy^jz^k>/<Xlm^2>  
c
c     and use expansion of Xlm to obtain
c
c     dinv(ijk,lm) = sum(i'j'k') di'j'k',lm <x^(i+i')y^(j+j')z^(k+k')> 
c     .              / <Xlm^2> 
c
c     <x^i y^j z^k> = 4Pi (i-1)!!(j-1)!!(k-1)!!/(i+j+k+1)!!
c
c     <Xlm^2> = 4Pi [(l+|m|)!/(l-|m|)!] / [(2-delta(m,0))(2l+1)]
c
c     (<> = integral over unit sphere)
c
      integer i, j, k, l, lp, m, mp, itri, ind, ijk, loff,
     $     ip, jp, kp, ijkp
      double precision numerator, denominator, sum
      double precision double_factorial, factorial
      external double_factorial, factorial
c
      itri(i,j)  = (i*(i-1))/2 + j
      ind(i,j,l) = itri(l-i+1,l-i-j+1) ! Index of x^i*y^j*z^(l-i-j)
      loff(l) = ((l*(l+1)*(l+2))/6)
c
      call xlm_coeff(lmax, d)
c
      call dfill((((lmax+1)*(lmax+2)*(lmax+3))/6)*(2*lmax+1)*(lmax+1),
     $     0.0d0, dinv, 1)
c
      do l = 0, lmax
         do i = l, 0, -1
            do j = l-i, 0, -1
               k = l-i-j
               ijk = loff(l) + ind(i,j,l)
               do lp = l, 0, -2
                  do ip = lp, 0, -1
                     do jp = lp-ip, 0, -1
                        kp = lp-ip-jp
                        if ( (mod(i+ip,2).eq.0) .and. 
     $                       (mod(j+jp,2).eq.0) .and. 
     $                       (mod(k+kp,2).eq.0) ) then
                        numerator = 
     $                       double_factorial(i+ip-1)*
     $                       double_factorial(j+jp-1)*
     $                       double_factorial(k+kp-1)/
     $                       double_factorial(i+ip+j+jp+k+kp+1)
                        ijkp = ind(ip,jp,lp)
                        do mp = 0, lp
                           denominator = factorial(lp+mp)/
     $                          (factorial(lp-mp)*(lp+lp+1))
*                           write(6,*) l, i, j, k
*                           write(6,*) lp, ip, jp, kp
*                           write(6,*) ijk, ijkp
*                           write(6,*) numerator, denominator, 
*     $                          d(ijkp, mp,lp),
*     $                          d(ijkp, -mp,lp)
                           if (mp .ne. 0) 
     $                          denominator = denominator*0.5d0
                           dinv(ijk,mp,lp) = dinv(ijk,mp,lp) + 
     $                          d(ijkp, mp,lp)*numerator/denominator
                           if (mp .ne. 0) 
     $                          dinv(ijk,-mp,lp) = dinv(ijk,-mp,lp) + 
     $                          d(ijkp,-mp,lp)*numerator/denominator
                        enddo
                        endif
                     enddo
                  enddo
               enddo
*               write(6,22) i, j, k, ijk
* 22            format(' Expansion for ', 3i5, 2x, i5)
*               do lp = l, 0, -2
*                  do mp = -lp,lp
*                     if (dinv(ijk,mp,lp) .ne. 0.0d0)
*     $                    write(6,1) mp, lp, l-lp, dinv(ijk,mp,lp)
* 1                   format(1x,2i5,3x,i5,3x,g20.10)
*                  end do
*               enddo
            enddo
         enddo
      enddo
c
c     Check what we can for paranoia.  The norm of a solid harmonic
c     is unchanged by converting to cartesians and back again.
c
      do l = 0, lmax
         do m = -l, l
            sum = 0.0d0
            do i = l, 0, -1
               do j = l-i, 0, -1
                  k = l-i-j
                  ijk = ind(i,j,l)
                  sum = sum + d(ijk,m,l)*dinv(loff(l)+ijk,m,l)
               enddo
            enddo
            if (abs(sum-1.0d0) .gt. 1d-10) stop 321
         enddo
      enddo
c
      end
      
      subroutine xlm_coeff(lmax, d)
      implicit none
      integer lmax
      double precision d(((lmax+1)*(lmax+2))/2, -lmax:lmax, 0:lmax)
c
c     compute the coefficients of cartesian polynomials for the
c     unnormalized real solid spherical harmonics up to order lmax
c
      integer i, j, k, l, m, itri, ind, lp1, ijk
      double precision twolp1

      itri(i,j)  = (i*(i-1))/2 + j
      ind(i,j,l) = itri(l-i+1,l-i-j+1) ! Index of x^i*y^j*z^(l-i-j)
c
      call dfill((((lmax+1)*(lmax+2))/2)*(2*lmax+1)*(lmax+1),
     $     0.0d0, d, 1)
c
      d(1, 0,0) = 1.0d0         ! Constant
c
      if (lmax .ge. 1) then
         d(1, 1,1) = 1.0d0      ! x
         d(2,-1,1) = 1.0d0      ! y
         d(3, 0,1) = 1.0d0      ! z
      end if
c
      do l = 1, lmax-1
         lp1 = l + 1
         twolp1 = 2*l + 1
c     
         do i = lp1, 0, -1
            do j = lp1-i, 0, -1
               k = lp1 - i - j
c
               ijk = ind(i,j,lp1)
c
*     q( lp1, lp1) = twolp1*(x*q( l,l) - y*q(-l,l))
               if (i.gt.0) d(ijk,lp1,lp1) =
     $              d(ijk,lp1,lp1) + twolp1*
     $              d(ind(i-1,j,l),l,l)
               if (j.gt.0) d(ijk,lp1,lp1) = 
     $              d(ijk,lp1,lp1) - twolp1*
     $              d(ind(i,j-1,l),-l,l)

*     q(-lp1, lp1) = twolp1*(x*q(-l,l) + y*q( l,l))
               if (i.gt.0) d(ijk,-lp1,lp1) =
     $              d(ijk,-lp1,lp1) + twolp1*
     $              d(ind(i-1,j,l),-l,l)
               if (j.gt.0) d(ijk,-lp1,lp1) = 
     $              d(ijk,-lp1,lp1) + twolp1*
     $              d(ind(i,j-1,l),l,l)

*     q(0,lp1) = (twolp1*z*q(0,l) - l*r2*q(0,l-1))/(l+1)
               if (k.gt.0) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) + twolp1*
     $              d(ind(i,j,l),0,l) / dble(l+1)
               if (i.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i-2,j,l-1),0,l-1) / dble(l+1)
               if (j.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i,j-2,l-1),0,l-1) / dble(l+1)
               if (k.gt.1) d(ijk,0,lp1) =
     $              d(ijk,0,lp1) - l*
     $              d(ind(i,j,l-1),0,l-1) / dble(l+1)

               do m = 1,l
*     q( m,lp1) = (twolp1*z*q( m,l) - (l+m)*r2*q( m,l-1))/(l-m+1)
*     q(-m,lp1) = (twolp1*z*q(-m,l) - (l+m)*r2*q(-m,l-1))/(l-m+1)

                  if (k.gt.0) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) +
     $                    twolp1*d(ind(i,j,l), m,l)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) +
     $                    twolp1*d(ind(i,j,l),-m,l)/dble(l-m+1)
                  end if
                  if (i.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i-2,j,l-1), m,l-1)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i-2,j,l-1),-m,l-1)/dble(l-m+1)
                  end if
                  if (j.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i,j-2,l-1), m,l-1)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i,j-2,l-1),-m,l-1)/dble(l-m+1)
                  end if
                  if (k.gt.1) then
                     d(ijk, m,lp1) = d(ijk, m,lp1) -
     $                    (l+m)*d(ind(i,j,l-1), m,l-1)/dble(l-m+1)
                     d(ijk,-m,lp1) = d(ijk,-m,lp1) -
     $                    (l+m)*d(ind(i,j,l-1),-m,l-1)/dble(l-m+1)
                  end if
               end do
c
c$$$               do m = -lp1,lp1
c$$$                  if (d(ijk,m,lp1) .ne. 0.0d0)
c$$$     $                 write(6,1) lp1, m, i, j, k, d(ijk,m,lp1)
c$$$ 1                format(1x,2i5,3x,3i3,3x,g20.10)
c$$$               end do
            end do
         end do
      end do
c$$$      stop 31
c
      end
      subroutine xlm_init
      implicit none
#include "xlm.fh"
c
c     initialize data for unnormalized solid spherical harmonics
c
c     b(m,l) = e(m) / (l + |m|)! , |m| <= l, 0 otherwise
c     c(m,l) = e(m) * (l - |m|)! , |m| <= l, 0 otherwise
c     twobc  = 2.0d0*b*c
c     rtwobc = 1.0d0/twobc
c     bcp0(m,l) = (b(m,l)*c(-m,l) + b(-m,l)*c(m,l))*phase(l+m) , m > 0
c     .           b(0,l)*c(0,l)*phase(l), m = 0
c     e(m) = -1 if m > 0 and m odd
c             1 otherwise
c     s(m) = -1 if m < 0
c             0 if m = 0
c             1 if m > 0
c     phase(m) = -1 if m odd  
c              =  1 if m even
c     qn(m,l) = normalization constant required by Condon & Shortley
c               so that the unnormalized xml satisfy
c               qn( m,l) * x( m,l) = r^l * Re (Yml) m >= 0
c               qn(-m,l) * x(-m,l) = r^l * Im (Yml) m >= 0
c               |m| <= l, 0 otherwise
c     rqn = 1.0d0 / qn
c     ee(m2,m1) = used to accellerate multipole_to_local
c     w, v are workspace
c     fastdiv(l) = 1.0d0/l, l=1,max2l+1
c
      integer l, m
      double precision pi4
      integer m1, m2, m1m2, m1mm2
      double precision pe1, pe1m

      numphi = 0                ! For caching rotation info
      numtheta = 0
      call ifill(maxhash, -1, iphis, 1)
      call ifill(maxhash, -1, ithetas, 1)

      do l = 1, maxl2+1
         fastdiv(l) = 1.0d0/dble(l)
      end do
      do l = 0, maxl
         do m = -maxl, maxl
            qn(m,l) = 0.0d0
            w(m,l) = 0.0d0
            v(m,l) = 0.0d0
         end do
      end do
      do l = 0, maxl2
         do m = -maxl2, maxl2
            b(m,l) = 0.0d0
            c(m,l) = 0.0d0
         end do
      end do
      do m = -maxl2, 0
         e(m) = 1.0d0
      end do
      do m = 1, maxl2, 2
         e(m) = -1.0d0
         if (m.lt.maxl2) e(m+1) = 1.0d0
      end do
      do m = -maxl2, maxl2
         s(m) = sign(1.0d0, dble(m))
      end do
      s(0) = 0.0d0
      phase(-maxl4) = 1.0d0
      do m = 1-maxl4, maxl4
         phase(m) = -phase(m-1)
      end do
c
      b(0,0) = 1.0d0
      c(0,0) = 1.0d0
      do l = 1, maxl2
         b(0,l) = b(0,l-1)/dble(l)
         c(0,l) = c(0,l-1)*dble(l)
         do m = 1, l
            b( m, l) = -b(m-1, l) / dble(l+m)
            b(-m, l) = abs(b(m,l))
            c(m, l) = -c(m-1, l) / dble(l-m+1)
            c(-m,l) = abs(c(m,l))
         end do
      end do
c
      do l = 0, maxl
         bcp0(0,l) = b(0,l)*c(0,l)*phase(l)
         do m = 1, l
            bcp0(m,l) = (b(m,l)*c(-m,l) + b(-m,l)*c(m,l))*phase(l+m)
         end do
      end do
c
      do l = 0, maxl
         do m = -l,l
            twobc(m,l) = 2.0d0*b(m,l)*c(m,l)
            rtwobc(m,l) = 1.0d0/twobc(m,l)
         enddo
      enddo
c
c     This defintion of the normalization constant is now correct
c     and agrees with numerical quadrature of the functions.
c     It also yields unitary rotation matrices.  
c
c     ANY changes must be checked for impact everywhere qn is used
c     since the correct normalization and unitary property 
c     are relied upon.
c
      pi4 = 16.0d0*atan(1.0d0)
      do l = 0, maxl
         do m = -l, l
            qn( m,l) = sqrt(2.0d0*(2*l+1)*b(m,l)*c(m,l)/pi4)*e(m)
            rqn(m,l) = 1.0d0/qn(m,l)
         end do
         qn(0,l) = qn(0,l) / sqrt(2.0d0)
         rqn(0,l) = 1.0d0/qn(0,l)
      end do
c
      do m1 = 0, maxl
         pe1 = phase(m1)*e(m1)
         pe1m= phase(-m1)
*         if (m1 .eq. 0) pe1 = pe1 * 0.5d0
*         if (m1 .eq. 0) pe1m = pe1m * 0.5d0
         do m2 = 0, maxl
            m1m2 = m1+m2
            m1mm2 = m1-m2
            ee( m2, m1) = (e(m2)*pe1*e(m1m2)+pe1m*e(-m1m2))*e(m1)
            ee(-m2, m1) = (pe1*e(m1mm2)+e(m2)*pe1m*e(-m1mm2))*e(m1)
         end do
      end do
c
*      write(6,*) ' b coefficients'
*      call xlm_print(maxl, b, maxl)
*      write(6,*) ' c coefficients'
*      call xlm_print(maxl, c, maxl)
*      write(6,*) ' qn coefficients'
*      call xlm_print(maxl, qn, maxl)
c
      end
      double precision function xlm_combin(n,m)
      implicit none
      integer n, m
c
c     Return n!/((n-m)!m!)
c
      double precision result
      integer i
c
      result = 1.0d0
      do i = n-m+1, n
         result = result * dble(i)
      enddo
      do i = 1, m
         result = result / dble(i)
      enddo
c
      xlm_combin = result
c
      end
      subroutine xlm_apply_z_y_rotation(lmax,alpha,beta,q,lq)
      implicit none
#include "xlm.fh"
      integer lmax,lq
      double precision alpha, beta, q(-lq:lq,0:lmax)
c
c     Apply a rotation of alpha radians about Z-axis followed
c     by a rotation of beta radians about the new Y-axis.
c
      call xlm_apply_z_rotation(lmax,alpha,q,lq)
      call xlm_apply_y_rotation(lmax,beta,q,lq)
c
      end
      subroutine xlm_apply_z_rotation(lmax,alpha,q,lq)
      implicit none
#include "xlm.fh"
      integer lmax,lq
      double precision alpha, q(-lq:lq,0:lmax)
c
c     Apply a rotation about the z axis of alpha radians
c     to the coefficients of unnormalized solid harmonics
c
      double precision cosa(0:maxl), sina(0:maxl), xp, xm
      integer l, m
c
      do m = 0,lmax
         cosa(m) = cos(m*alpha)
         sina(m) = sin(m*alpha)
      enddo
c
      do l = 1, lmax
         do m = 1, l
            xp = cosa(m)*q(m,l) + sina(m)*q(-m,l)
            xm =-sina(m)*q(m,l) + cosa(m)*q(-m,l)
            q( m,l) = xp
            q(-m,l) = xm
         enddo
      enddo
c
      end
      subroutine xlm_apply_y_rotation(lmax,beta,q,lq)
      implicit none
#include "xlm.fh"
      integer lmax,lq
      double precision beta, q(-lq:lq,0:lmax)
c
c     Apply a rotation about the y axis of beta radians
c     to the coefficients of unnormalized solid harmonics
c
      double precision dp((maxl+1)*(maxl+1)), 
     $     dm(maxl*maxl), qq(-maxl:maxl)
      integer l, m, mp, ind
c
      do l = 0, lmax
         call xlm_y_rotation_matrix(l,beta,dp,dm)
         do m = -l, l
            qq(m) = q(m,l)*rqn(m,l) ! Precondition
            q(m,l) = 0.0d0
         enddo
         ind = 1
         do mp = 0, l
            do m = 0, l
               q(m,l) = q(m,l) + dp(ind)*qq(mp)
               ind = ind + 1
            enddo
         enddo
         if (l .gt. 0) then
            ind = 1
            do mp = -l,-1
               do m = -l,-1
                  q(m,l) = q(m,l) + dm(ind)*qq(mp)
                  ind = ind + 1
               enddo
            enddo
         endif
         do m = -l, l
            q(m,l) = q(m,l)*qn(m,l) ! Undo preconditioning
         enddo
      enddo
c
      end
      subroutine xlm_y_rotation_matrix(l,beta,dp,dm)
      implicit none
#include "xlm.fh"
c
      integer l
      double precision beta, dp(0:l,0:l),dm(l,l)
c
c     Return in dp and dm the matrices corresponding to the
c     operator rotating about the y axis by angle beta.
c
c     For numerical stability, the rotation matrices are
c     PRECONDITIONED with qn(m,l), the normalization constant
c     for the funciton.
c
c     Note that the preconditioned matrices are unitary.
c
c     Thus,
c
c     qn(m,l)*Xlm (theta+angle) = sum(m') Xlm' (theta)*qn(m',l)*dp(m',m)
c     qn(-m,l)*Xl-m(theta+angle) = sum(m') Xl-m'(theta)*qn(-m',l)*dp(-m',-m)
c
c     Or, to apply the operator to a set of coefficients consider
c
c     f(theta) = sum(m) Cm Xlm(theta) 
c     .        = sum(m) (Cm/qn(m))*(Xlm(theta)*qn(m))
c     .        = sum(m,m') (Cm/qn(m))*(Xlm'(theta+angle)*qn(m')*Dm'm)
c     .        = sum(m') Xlm'(theta+angle)*qn(m')*[sum(m) Dm'm*(Cm/qn(m))]
c     
c     i.e., the coeff of Xlm'(theta+beta) is qn(m')*sum(m) Dm'm*(Cm/qn(m))
c     which are formed by first dividing the coeffs by qn(m), multiplying 
c     by Dm'm and then multiplying by qn(m').
c
c     This is much more numerically stable than not scaling.
c
c     This routine is meant to be correct and accurate, but there
c     are MUCH faster ways of generating these matrices.
c     
      integer k, info, ipiv(maxl+1), m
      double precision theta, pi, x, y, z, rx, ry, rz, phi
      double precision a(0:maxl,0:maxl)
c
      pi = 4.0d0*atan(1.0d0)
c
      do k = 0, l
         theta = dble(k)*pi/dble(l+1)
         z = cos(theta+beta)
         x = sin(theta+beta)
         y = 0.0d0
         call xlm(l, x, y, z, w, maxl)
         do m = 0, l
            dp(k,m) = w(m,l)*qn(m,l) ! Condition with qn
         enddo
         z = cos(theta)
         x = sin(theta)
         y = 0.0d0
         call xlm(l, x, y, z, w, maxl)
         do m = 0, l
            a(k,m) = w(m,l)*qn(m,l) ! Condition with qn
         enddo
      enddo
c
      call dgesv(l+1,l+1,a,maxl+1,ipiv,dp,l+1,info)
      if (info .ne. 0) call errquit('xlm_y_rot_mat: Pinfo ', info)
c
c$$$      do m = 0, l
c$$$         dp(m,m) = dp(m,m) - 1.0d0
c$$$      enddo
c$$$      write(6,*) ' DP '
c$$$      call doutput(dp,1,l+1,1,l+1,l+1,l+1,1)
c$$$      call dfill((l+1)**2,0.0d0,dp,1)
c$$$      do m = 0, l
c$$$         dp(m,m) = dp(m,m) + 1.0d0
c$$$      enddo
c
      if (l .gt. 0) then
         phi = pi/dble(2*l+1)     ! -ve m components are zero at phi=0
         do k = 1, l
            theta = dble(k)*pi/dble(l+1)
            z = cos(theta)
            x = sin(theta)*cos(phi)
            y = sin(theta)*sin(phi)
            call xlm(l, x, y, z, w, maxl)
            do m = 1, l
               a(k,l-m+1) = w(-m,l)*qn(-m,l)
            enddo
c     
            rx = x*cos(-beta) - z*sin(-beta)
            ry = y
            rz = z*cos(-beta) + x*sin(-beta)
            call xlm(l, rx, ry, rz, w, maxl)
            do m = 1, l
               dm(k,l-m+1) = w(-m,l)*qn(-m,l)
            enddo
         enddo
c
         call dgesv(l,l,a(1,1),maxl+1,ipiv,dm,l,info)
         if (info .ne. 0) call errquit('xlm_y_rot_mat: Minfo ', info)
c
c$$$         do m = 1, l
c$$$            dm(m,m) = dm(m,m) - 1.0d0
c$$$         enddo
c$$$         write(6,*) ' DM '
c$$$         call doutput(dm,1,l,1,l,l,l,1)
c$$$         call dfill(l*l,0.0d0,dm,1)
c$$$         do m = 1, l
c$$$            dm(m,m) = dm(m,m) + 1.0d0
c$$$         enddo
      endif
c
      end
      
      subroutine xlm_translate(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      double precision x, y, z
      integer lmax, lq
      double precision q(-lq:lq,0:lmax)
      integer lbig, l, m, mmp, mp, mmpa, lp, llp, mpa
      double precision sump, summ, bfac
c
c     given a set of unnormalized multipoles about the origin in q
c     translate the center to (x,y,z).
c
      if (lmax .gt. maxl) then
         write(6,*) ' parameter exceeded', lmax, maxl
         stop 1
      end if
c
      do l = lmax, 0, -1
         do m = -l, l
            if (abs(q(m,l)) .gt. 0.0d0) goto 10
         enddo
      enddo

 10   lbig = l
      if (lbig .lt. 0) then
         do l = 0, lmax
            do m = -l, l
               q(m,l) = 0.0d0
            enddo
         enddo
         return
      endif
c
      call xlm(lmax, -x, -y, -z, w, maxl)
c
c     transfer the charge, then everything else
c
      do l = 0, lmax
         do m = -l, l
*            q(m,l) = q(m,l) ! UH ?
*            w(m,l) = w(m,l)
            v(m,l) = w(m,l) * q(0,0) * b(abs(m),l) * b(0,0)
         end do
      end do
c
      do l = 1, lmax
         do m = 0, l
            sump = 0.0d0
            summ = 0.0d0
            do lp = 1, min(l,lbig)
               llp = l - lp
               do mp = -lp, lp
                  mmp = m - mp
                  mmpa = abs(mmp)
                  if (mmpa .le. llp) then
                     bfac = b(mmp,llp)*b(mp,lp)
                     mpa = abs(mp)
                     sump = sump + (w(mmpa,llp)*q(mpa,lp) - 
     $                    w(-mmpa,llp)*q(-mpa,lp)*s(mmp)*s(mp))*bfac
                     summ = summ + (w(mmpa,llp)*q(-mpa,lp)*s(mp) +
     $                    w(-mmpa,llp)*q(mpa,lp)*s(mmp))*bfac
                  end if
               end do
            end do
            v(-m,l) = v(-m,l) + summ
            v( m,l) = v( m,l) + sump
         end do
      end do
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = v(m,l) / b(abs(m),l)
         end do
      end do
c
      end
      subroutine xlm_local_translate(lmax, x, y, z, q, lq)
      implicit none
#include "xlm.fh"
c
      integer lmax, lq
      double precision x, y, z, q(-lq:lq,0:lmax)
c
c     Given a local spherical taylor series for the potential 
c     translate the center of the expansion to (x,y,z).  
c
c     If the routine xlm_translate does
c     .   Qml(r-a) <= sum(m'l')Tmlm'l'(-a)Qm'l'(r)
c     this routine does
c     .   Qml(r-a) <= sum(m'l')Qm'l'(r)Tm'l'ml(a)
c
      integer l, m, mp, mmp, mmpa, mpa, lp, llp
      double precision bfac
c
      if (lmax .gt. maxl) then
         write(6,*) ' parameter exceeded', lmax, maxl
         stop 1
      end if
c
      call xlm(lmax, x, y, z, w, maxl)
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = q(m,l) / b(abs(m),l)
            v(m,l) = 0.0d0
         end do
      end do
c
      do l = 0, lmax
         do m = 0, l
            do lp = 0, l
               llp = l - lp
               do mp = -lp, lp
                  mmp = m - mp
                  mmpa = abs(mmp)
                  if (mmpa .le. llp) then
                     bfac = b(mmp,llp)*b(mp,lp)
                     mpa = abs(mp)

*                     Tm,l,mpa,lp  +=  w(mmpa,llp)*bfac
*                     Tm,l,-mpa,lp += -w(-mmpa,llp)*s(mmp)*s(mp)*bfac
*                     T-m,l,-mpa,lp += w(mmpa,llp)*s(mp)*bfac
*                     T-m,l,mpa,lp  += w(-mmpa,llp)*s(mmp))*bfac

                     v( mpa,lp) = v( mpa,lp) + bfac*(
     $                    w(mmpa,llp)*q(m,l) + 
     $                    w(-mmpa,llp)*s(mmp)*q(-m,l))
                     v(-mpa,lp) = v(-mpa,lp) + bfac*(
     $                    w(mmpa,llp)*s(mp)*q(-m,l) - 
     $                    w(-mmpa,llp)*s(mmp)*s(mp)*q(m,l))

                  end if
               end do
            end do
         end do
      end do
c
      do l = 0, lmax
         do m = -l, l
            q(m,l) = v(m,l)
         end do
      end do
c
      end
      double precision function xlm_coulomb(lmax, ax, ay, az, aq, la,
     $     bx, by, bz, bq, lb)
      implicit none
#include "xlm.fh"
c
      integer lmax, la, lb
      double precision ax, ay, az, bx, by, bz
      double precision aq(-la:la, 0:lmax), bq(-lb:lb, 0:lmax)
c
c     compute the coulomb interaction between two unnormalized
c     real solid multipole expansions at centers a and b
c
      double precision xx, yy, zz, rsq, r1, r12, energy, ei, ereal
      integer l1, l2, m1, m2, l2max, m1a, m2a, l12, m12, m12a
      double precision bcpr, bq_plus, bq_minus, aq_plus, aq_minus
      double precision qnorm, eimag
c
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(2*lmax, xx, yy, zz, u, maxl2)
c
      energy = 0.0d0
      ei = 0.0d0
      r1 = dsqrt(rsq)
      do l1 = 0, lmax
         l2max = lmax        ! Full summation
         do m1 = -l1, l1
            qnorm = abs(aq(m1,l1)) + abs(aq(-m1,l1))
            if (qnorm .ne. 0.0d0) then
               m1a = abs(m1)
               ereal = 0.0d0
               eimag = 0.0d0
               r12 = r1 * b(m1,l1) ! Included b here to keep no.s small
               do l2 = 0, l2max
                  do m2 = -l2, l2
                     m2a = abs(m2)
                     l12 = l1+l2
                     m12 = -(m1+m2)
                     m12a = abs(m12)
                     bcpr = (b(m2,l2)*r12)*(c(m12,l12)*phase(l2+m12))
                     bq_plus = (bq(m2a,l2)*bcpr)
                     bq_minus= (bq(-m2a,l2)*bcpr)
                     ereal = ereal + (bq_plus*u(m12a,l12) -
     $                    bq_minus*u(-m12a,l12)*s(m2)*s(m12))
                     eimag = eimag + (bq_minus*s(m2)*u(m12a,l12)+
     $                    bq_plus*u(-m12a,l12)*s(m12))
                  end do
                  r12 = r12 * rsq
               end do
               aq_plus = (aq(m1a,l1))
               aq_minus= (aq(-m1a,l1))
               energy = energy + 
     $              (ereal*aq_plus - eimag*aq_minus*s(m1))
               ei = ei + (ereal*aq_minus*s(m1)+eimag*aq_plus)
            end if
         end do
         r1 = r1 * rsq
      end do
*      if (abs(ei) .gt. 1.0d-10) then
*         write(6,*) ' xlmcoul: imaginary part of potential is large',
*     $        energy, ei
*      end if
c
      xlm_coulomb = energy
c
      end
      subroutine xlm_print(lmax, q, lq)
      implicit double precision (a-h, o-z)
      dimension q(-lq:lq, 0:lmax)
      parameter(mchunk = 6)
c
      do l = 0, lmax
         sum = 0.0d0
         do m = -l,l
            sum = sum + abs(q(m,l))
         end do
         if (sum .le. 1e-40) then
            write(6,3) l
 3          format(' l = ',i2,' is zero')
         else
            write(6,1) l
 1          format(' l = ',i2)
            do mlow = -l, l, mchunk
               mhi = min(l, mlow+mchunk-1)
               write(6,2) (q(m,l),m=mlow,mhi)
 2             format('      ',1p,6d12.4)
            end do
         end if
      end do
c
      end
      double precision function xlm_potential(ax, ay, az, 
     $     bx, by, bz, bq, lb, lmax)
      implicit none
#include "xlm.fh"
c
      integer lb, lmax
      double precision ax, ay, az, bx, by, bz
      double precision bq(-lb:lb, 0:lmax)
c
c     compute the coulomb potential at A from the unnormalized 
c     real-solid multipole distribution at B.
c
      double precision xx, yy, zz, rsq
      double precision energy, r12
      integer l, m
c
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax, xx, yy, zz, u, maxl2)
c
      energy = 0.0d0
      r12 = dsqrt(rsq)
      do l = 0, lmax
         energy = energy + bcp0(0,l)*r12*bq(0,l)*u(0,l)
         do m = 1, l
            energy = energy + (bq(m,l)*u(m,l) +
     $           bq(-m,l)*u(-m,l))*bcp0(m,l)*r12
         end do
         r12 = r12 * rsq
      end do
c
      xlm_potential = energy
c
      end
      double precision function xlm_potential_test(ax, ay, az, 
     $     bx, by, bz, bq, lb, lmax)
      implicit none
#include "xlm.fh"
c
      integer lb, lmax
      double precision ax, ay, az, bx, by, bz
      double precision bq(-lb:lb, 0:lmax)
c
c     This version of xlm_potential just to make sure
c     the formula in the paper is correct!  It is.
c
      double precision xx, yy, zz, rsq
      double precision energy, r12, fourpi, sum
      integer l, m
c
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax, xx, yy, zz, u, maxl2)
      fourpi = 16.0d0*atan(1.0d0)
c
      energy = 0.0d0
      r12 = dsqrt(rsq)
      do l = 0, lmax
         sum = 0.0d0
         do m = -l, l
            sum = sum + bq(m,l)*qn(m,l)*qn(m,l)*u(m,l)
         end do
         energy = energy + sum*fourpi*r12/dble(l+l+1)
         r12 = -r12 * rsq
      end do
c
      xlm_potential_test = energy
c
      end
      double precision function xlm_local_potential(x, y, z, 
     $     aq, la, lmax)
      implicit none
#include "xlm.fh"
      double precision x, y, z
      integer la, lmax
      double precision aq(-la:la,0:lmax)
c
c     Given in aq() an expansion of the local potential compute
c     it at the given point.  Note, that if the potential is 
c     expanded about some point A, then the coordinates relative
c     to A (r-a) should be passed in.
c
      integer l, m
      double precision sum
c
      call xlm(lmax, x, y, z, w, maxl)
      sum = 0.0d0
      do l = 0, lmax
         do m = -l,l
            sum = sum + w(m,l)*aq(m,l)
         end do
      end do
c
      xlm_local_potential = sum
c
      end
      subroutine xlm_unpack(lmax, p, q, lq)
      implicit none
      integer lmax, lq
      double precision p(*), q(-lq:lq, 0:lmax)
c
c     Unpack the multipoles in p (which are stored contiguously
c     with no space between) into the easier to use array q()
c
      integer l, m, ind
c
      ind = 1
      do l = 0, lmax
         do m = -l, l
            q(m,l) = p(ind)
            ind = ind + 1
         end do
      end do
c
      end
      subroutine xlm_pack(lmax, p, q, lq)
      implicit none
      integer lmax, lq
      double precision p(*), q(-lq:lq, 0:lmax)
c
c     Pack the multipoles into p (which are stored contiguously
c     with no space between) from the easier to use array q()
c
      integer l, m, ind
c
      ind = 1
      do l = 0, lmax
         do m = -l, l
            p(ind) = q(m,l)
            ind = ind + 1
         end do
      end do
c
      end
      subroutine xlm_accumulate_to_packed(lmax, p, q, lq)
      implicit none
      integer lmax, lq
      double precision p(*), q(-lq:lq, 0:lq)
c
c     Add into the packed multipoles p (which are stored contiguously
c     with no space between) from the easier to use array q()
c
      integer l, m, ind
c
      ind = 1
      do l = 0, lmax
         do m = -l, l
            p(ind) = p(ind) + q(m,l)
            ind = ind + 1
         end do
      end do
c
      end
      subroutine xlm_multipole_to_local(lmax, ax, ay, az, aq, la,
     $     bx, by, bz, bq, lb)
      implicit none
#include "xlm.fh"
c     
      integer lmax, la, lb
      double precision ax, ay, az, bx, by, bz
      double precision aq(-la:la, 0:la), bq(-lb:lb, 0:lb)
c     
c     Put into aq() the multipolar potential from the unnormalized
c     real solid multipole expansion at center B computed at center A.
c     
c     I.e., the potential at a point r near A due to the multipoles
c     at B may be computed as  
c     
c     .     sum(lm) Xlm(r-a)AQlm
c
c     
      double precision xx, yy, zz, rsq, r12, ereal, eimag
      integer l1, l2, m1, m2, l2max, l12, m12, m12a
      integer m1m2, m1mm2
      double precision urpp, urmm, uipp, uimm
**      integer flops
      double precision ww(0:maxl,0:maxl)
      double precision wv(0:maxl,0:maxl)
      double precision vw(0:maxl,0:maxl)
      double precision vv(0:maxl,0:maxl)
c
      do l1 = 0, la
         do m1 = -l1,l1
            aq(m1,l1) = 0.0d0
         enddo
      enddo
c
      do l2 = lmax,0,-1
         do m2 = -l2,l2
            if (abs(bq(m2,l2)) .gt. 0.0d0) goto 10
         enddo
      enddo
 10   l2max = l2
c
      if (l2max .lt. 0) return
c         
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax+l2max, xx, yy, zz, u, maxl2)
c
**      flops = 0
c
      do l2 = 0, l2max
         do m2 = 0, l2
            w(m2,l2) = bq( m2,l2)*b(m2,l2)*phase(l2+m2)*e(m2)
            v(m2,l2) = bq(-m2,l2)*b(m2,l2)*phase(l2+m2)*s(m2)*e(m2)
         enddo
         w(0,l2) = w(0,l2)*0.5d0
         v(0,l2) = v(0,l2)*0.5d0
      enddo
c
      do l12 = 0, lmax+l2max
         do m12 = -l12,l12
            m12a = abs(m12)
            ur(m12,l12) = u( m12a,l12)*c(-m12,l12)*e(m12)
            ui(m12,l12) = u(-m12a,l12)*c(-m12,l12)*s(-m12)*e(m12)
         end do
      end do
c
      r12 = sqrt(rsq)
      do l12 = 0, lmax+l2max
*         call dfill((maxl+1)**2,0.0d0,ww,1)
*         call dfill((maxl+1)**2,0.0d0,wv,1)
*         call dfill((maxl+1)**2,0.0d0,vw,1)
*         call dfill((maxl+1)**2,0.0d0,vv,1)
         do m1 = 0, min(l12,lmax)
            do m2 = 0, min(m1,l12-m1,l2max)
               m1m2 = m1+m2
               m1mm2 = m1-m2
               urpp = ur(m1m2, l12)*ee( m2,m1)
               urmm = ur(m1mm2,l12)*ee(-m2,m1)
               ww(m2,m1) = urpp + urmm
               ww(m1,m2) = ww(m2,m1)
               vv(m2,m1) = urpp - urmm
               vv(m1,m2) = vv(m2,m1)
c
               uipp = ui(m1m2, l12)*ee( m2,m1)
               uimm = ui(m1mm2,l12)*ee(-m2,m1)
               wv(m2,m1) = uimm - uipp
               vw(m1,m2) = -wv(m2,m1)
               vw(m2,m1) = uimm + uipp
               wv(m1,m2) = -vw(m2,m1)
            enddo
**            flops = flops + 8*(1+min(m1,l12-m1,l2max))
         enddo
c$$$         lll = min(l12+1,lmax+1)
c$$$         write(6,*) l12, ' ww '
c$$$         call doutput(ww, 1, lll, 1, lll, maxl+1, maxl+1, 1)
c$$$         write(6,*) l12, ' vv '
c$$$         call doutput(vv, 1, lll, 1, lll, maxl+1, maxl+1, 1)
c$$$         write(6,*) l12, ' wv '
c$$$         call doutput(wv, 1, lll, 1, lll, maxl+1, maxl+1, 1)
c$$$         write(6,*) l12, ' vw '
c$$$         call doutput(vw, 1, lll, 1, lll, maxl+1, maxl+1, 1)
         do l1 = max(0,l12-l2max), min(l12,lmax)
            l2 = l12 - l1
            ereal = ww(0,0)*w(0,l2)
            do m2 = 1, l2
               ereal = ereal+ww(m2,0)*w(m2,l2)+wv(m2,0)*v(m2,l2)
            end do
            aq( 0,l1) = aq( 0,l1) + ereal*r12
**            flops = flops + 3 + l2*4
            do m1 = 1, l1
               ereal = ww(0,m1)*w(0,l2)
               eimag = vw(0,m1)*w(0,l2)
               do m2 = 1, l2
                  ereal = ereal+ww(m2,m1)*w(m2,l2)+wv(m2,m1)*v(m2,l2)
                  eimag = eimag+vw(m2,m1)*w(m2,l2)+vv(m2,m1)*v(m2,l2)
               end do
               aq( m1,l1) = aq( m1,l1) + ereal*r12
               aq(-m1,l1) = aq(-m1,l1) - eimag*r12
            end do
**            flops = flops + l1*(6+l2*8)
         end do
         r12 = r12 * rsq
      end do
c
      do l1 = 0,lmax
         aq( 0,l1) = aq(0,l1)*b(0,l1)*0.5d0
         do m1 = 1, l1
            aq( m1,l1) = aq( m1,l1)*b(m1,l1)
            aq(-m1,l1) = aq(-m1,l1)*b(m1,l1)*s(m1)
         enddo
      enddo
c
c
**       write(6,*) ' oldflops ', flops
*       stop
c
      end
      subroutine xlm_multipole_to_local_z(lmax, az, aq, la, bq, lb)
      implicit none
#include "xlm.fh"
c     
      integer lmax, la, lb
      double precision az
      double precision aq(-la:la, 0:la), bq(-lb:lb, 0:lb)
c     
c     Put into aq() the multipolar potential from the unnormalized
c     real solid multipole expansion at center B computed at center A.
c
c     RESTRTICTED TO Z ONLY TRANSLATION
c     
c     I.e., the potential at a point r near A due to the multipoles
c     at B may be computed as  
c     
c     .     sum(lm) Xlm(r-a)AQlm
c     
      double precision r, r12, r12ww
      integer l1, l2, m1, m2, l2max, l12, minl1l2
**      integer flops
      double precision ww, rsign
c
      do l1 = 0, la
         do m1 = -l1,l1
            aq(m1,l1) = 0.0d0
         enddo
      enddo
c
      do l2 = lmax,0,-1
         do m2 = -l2,l2
            if (abs(bq(m2,l2)) .gt. 0.0d0) goto 10
         enddo
      enddo
 10   l2max = l2
***      l2max = lmax
c
      if (l2max .lt. 0) return
c
      do l2 = 0, l2max
         do m2 = -l2, l2
            w(m2,l2) = bq( m2,l2)*b(m2,l2)*phase(l2+m2)
         enddo
**         flops = flops + 2*(2*l2+1)
      enddo
c
      r = 1.0d0/az
      r12 = r
c
c     l12 = 0 
c
      ww = 2.0d0
      aq(0,0) = aq(0,0) + r12*ww*w(0,0)
      r12 = r12 * r
**      flops = flops + 4
c
c     l12 = 1
c
      if (lmax+l2max .ge. 1) then
         ww = -ww
         r12ww = r12*ww
         aq(0,0) = aq(0,0) + r12ww*w(0,1) ! l1=0, l2=1
         aq(0,1) = aq(0,1) + r12ww*w(0,0) ! l1=1, l2=0
         r12 = r12 * r
**         flops = flops + 10
      endif
c
      do l12 = 2, lmax+l2max
         ww = -ww * dble(l12) ! (-1)^l * 2 * l!
         do l1 = max(0,l12-l2max), min(l12,lmax)
            l2 = l12 - l1
            minl1l2 = min(l1,l2)
            r12ww = r12*ww
            do m1 = -minl1l2, minl1l2
               aq( m1,l1) = aq( m1,l1) + r12ww*w(m1,l2)
            end do
**            flops = flops + 2 + 2*(2*minl1l2 + 1)
         end do
         r12 = r12 * r
**         flops = flops + 3
      end do
c
      rsign = sign(1.0d0,r)
      do l1 = 0,lmax
         aq( 0,l1) = aq(0,l1)*b(0,l1)*0.5d0*rsign
         do m1 = 1, l1
            aq( m1,l1) = aq( m1,l1)*b(m1,l1)*rsign
            aq(-m1,l1) = aq(-m1,l1)*b(m1,l1)*rsign*phase(m1)
         enddo
**         flops = flops + 3*l1 + 2
      enddo
c
c
**      write(6,*) ' flops ', lmax, flops
c
      end
      subroutine xlm_multipole_to_local_swap(lmax,
     $     ax, ay, az, aq, ap, la,
     $     bx, by, bz, bq, bp, lb)
      implicit none
#include "xlm.fh"
c     
      integer lmax, la, lb
      double precision ax, ay, az, bx, by, bz
      double precision aq(-la:la, 0:la), bq(-lb:lb, 0:lb)
      double precision ap(-la:la, 0:la), bp(-lb:lb, 0:lb)
c     
c     Put into ap() the multipolar potential from the unnormalized
c     real solid multipole expansion at center B computed at center A
c     and similarly pu the potential from aq() into bp()
c     
c     I.e., the potential at a point r near A due to the multipoles
c     at B may be computed as  
c     
c     .     sum(lm) Xlm(r-a)APlm
c
c     and vice versa
c     
      double precision xx, yy, zz, rsq, r12, r12ph
      double precision ereala, eimaga, erealb, eimagb 
      integer l1, l2, m1, m2, l12, m12, m12a
      integer m1m2, m1mm2, lmax2
      double precision urpp, urmm, uipp, uimm, ce, bpe
      integer flops
      double precision ww(0:maxl,0:maxl)
      double precision wv(0:maxl,0:maxl)
      double precision vw(0:maxl,0:maxl)
      double precision vv(0:maxl,0:maxl)
      double precision wa(0:maxl,0:maxl), va(0:maxl,0:maxl)
      double precision wb(0:maxl,0:maxl), vb(0:maxl,0:maxl)
c
      lmax2 = lmax+lmax
c
      do l1 = 0, la
         do m1 = -l1,l1
            ap(m1,l1) = 0.0d0
         enddo
      enddo
      do l1 = 0, lb
         do m1 = -l1,l1
            bp(m1,l1) = 0.0d0
         enddo
      enddo
c         
      xx = bx - ax
      yy = by - ay
      zz = bz - az
      rsq = 1.0d0 / (xx*xx + yy*yy + zz*zz)
      call xlm(lmax2, xx, yy, zz, u, maxl2)
c
      flops = 0
c
      do l2 = 0, lmax
         do m2 = 0, l2
            bpe = b(m2,l2)*phase(l2+m2)*e(m2)
            wa(m2,l2) = aq( m2,l2)*bpe
            va(m2,l2) = aq(-m2,l2)*(bpe*s(m2))
            wb(m2,l2) = bq( m2,l2)*bpe
            vb(m2,l2) = bq(-m2,l2)*(bpe*s(m2))
         enddo
         wa(0,l2) = wa(0,l2)*0.5d0
         wb(0,l2) = wb(0,l2)*0.5d0
      enddo
c
      do l12 = 0, lmax2
         do m12 = -l12,l12
            m12a = abs(m12)
            ce = c(-m12,l12)*e(m12)
            ur(m12,l12) = u( m12a,l12)*ce
            ui(m12,l12) = u(-m12a,l12)*ce*s(-m12)
         end do
      end do
c
      r12   = sqrt(rsq)
      r12ph = r12
      do l12 = 0, lmax2
         do m1 = 0, min(l12,lmax)
            do m2 = 0, min(m1,l12-m1,lmax)
               m1m2 = m1+m2
               m1mm2 = m1-m2
               urpp = ur(m1m2, l12)*ee( m2,m1)
               urmm = ur(m1mm2,l12)*ee(-m2,m1)
               ww(m2,m1) = urpp + urmm
               ww(m1,m2) = ww(m2,m1)
               vv(m2,m1) = urpp - urmm
               vv(m1,m2) = vv(m2,m1)
c
               uipp = ui(m1m2, l12)*ee( m2,m1)
               uimm = ui(m1mm2,l12)*ee(-m2,m1)
               wv(m2,m1) = uimm - uipp
               vw(m1,m2) = -wv(m2,m1)
               vw(m2,m1) = uimm + uipp
               wv(m1,m2) = -vw(m2,m1)
            enddo
*            flops = flops + 8*(1+min(m1,l12-m1,lmax))
         enddo
         do l1 = max(0,l12-lmax), min(l12,lmax)
            l2 = l12 - l1
            ereala = ww(0,0)*wa(0,l2)
            erealb = ww(0,0)*wb(0,l2)
            do m2 = 1, l2
               ereala = ereala + ww(m2,0)*wa(m2,l2)+wv(m2,0)*va(m2,l2)
               erealb = erealb + ww(m2,0)*wb(m2,l2)+wv(m2,0)*vb(m2,l2)
            end do
            ap( 0,l1) = ap( 0,l1) + erealb*r12
            bp( 0,l1) = bp( 0,l1) + ereala*r12ph
*            flops = flops + 3 + l2*4
            do m1 = 1, l1
               ereala = ww(0,m1)*wa(0,l2)
               eimaga = vw(0,m1)*wa(0,l2)
               do m2 = 1, l2
                  eimaga =eimaga+vw(m2,m1)*wa(m2,l2)+vv(m2,m1)*va(m2,l2)
                  ereala =ereala+ww(m2,m1)*wa(m2,l2)+wv(m2,m1)*va(m2,l2)
               enddo
               bp( m1,l1) = bp( m1,l1) + ereala*r12ph
               bp(-m1,l1) = bp(-m1,l1) - eimaga*r12ph
               erealb = ww(0,m1)*wb(0,l2)
               eimagb = vw(0,m1)*wb(0,l2)
               do m2 = 1, l2    ! SUN slower with merged loops ?
                  erealb =erealb+ww(m2,m1)*wb(m2,l2)+wv(m2,m1)*vb(m2,l2)
                  eimagb =eimagb+vw(m2,m1)*wb(m2,l2)+vv(m2,m1)*vb(m2,l2)
               end do
               ap( m1,l1) = ap( m1,l1) + erealb*r12
               ap(-m1,l1) = ap(-m1,l1) - eimagb*r12
            end do
*            flops = flops + l1*(6+l2*8)
         end do
         r12 = r12 * rsq
         r12ph = -r12ph * rsq
      end do
c
      do l1 = 0,lmax
         ap( 0,l1) = ap( 0,l1)*b(0,l1)*0.5d0
         bp( 0,l1) = bp( 0,l1)*b(0,l1)*0.5d0
         do m1 = 1, l1
            ap( m1,l1) = ap( m1,l1)*b(m1,l1)
            ap(-m1,l1) = ap(-m1,l1)*b(m1,l1)*s(m1)
            bp( m1,l1) = bp( m1,l1)*b(m1,l1)
            bp(-m1,l1) = bp(-m1,l1)*b(m1,l1)*s(m1)
         enddo
      enddo
c
c
*       write(6,*) ' flops ', flops
*       stop
c
      end
      subroutine xlm_new_multipole_to_local(lmax, 
     $     xp, yp, zp, p, lp, xq, yq, zq, q, lq)
      implicit none
#include "xlm.fh"
c     
      integer lmax, lq, lp
      double precision xp, yp, zp, xq, yq, zq
      double precision q(-lq:lq,0:lmax), p(-lp:lp,0:lp)
c
      double precision theta, phi, r, x, y, z, zz, phisign, qqp, qqm
      double precision pi, thetasign
      double precision qq(-maxl:maxl,0:maxl), tmp(-maxl:maxl)
      integer l, m, mp, ihash, iphi, itheta, indp, indm
c
c     From the multipole-potential q() at (xq,yq,zq) compute
c     the local Taylor series potential at (xp,yp,zp), in p() 
c     using a rotation-based algorithm.  
c
c     Under the assumption that only a limited no. of angles
c     will be used, cache rotation matrix info.
c
      x = xp - xq
      y = yp - yq
      z = zp - zq
      r = sqrt(x*x + y*y + z*z)
c
c     Determine the required rotations
c
      pi = 4.0d0*atan(1.0d0)
      if (x.eq.0.0d0 .and. y.eq.0.0d0) then
         phi = 0.0d0
         theta = 0.0d0
         if (z .lt. 0.0d0) theta = pi
      else 
         zz = z/r
         if (abs(zz) .gt. 1.0d0) zz = sign(1.0d0,zz)
         theta = acos(zz)
         phi = atan2(y,x)
      endif
c      
c     Rotate the source multipoles.  The actual potential
c     is roughly 2*qml*bml*cml*xml/r^(2l+1), so before we
c     rotate we must include these additional factors.
c
c     Hash into the list of cached angles to get the rotation info.
c     The cached phis are all positive so need to track the sign.
c     The cached thetas are all 0-pi/2, so if phi is pi/2-pi then
c     must rotate with theta-pi and use negative translation.
c
      ihash = int(abs(phi*1000.0d0))
      if (ihash .gt. maxhash) call errquit('xlm_n_m_l: phi hash',ihash)
      if (iphis(ihash) .gt. 0) then
         if (abs(phis(ihash)-abs(phi)) .gt. 1d-12) then
            write(6,*) ' ihash phi phis ', ihash, phi, phis(ihash)
            call errquit('xlm_n_m_l: duplicate phi hash',0)
         endif
      else
         call xlm_cache_phi(phi)
      endif
      iphi = iphis(ihash)
      phisign = sign(1.0d0,phi)
c
      if (theta .le. 0.5d0*pi) then
         thetasign = 1.0d0
      else
         theta = theta - pi
         thetasign = -1.0d0
      endif
      ihash = int(abs(theta*1000.0d0))
      if (ihash .gt. maxhash)call errquit('xlm_n_m_l: theta hash',ihash)
      if (ithetas(ihash) .gt. 0) then
         if (abs(thetas(ihash)-abs(theta)) .gt. 1d-12) then
            write(6,*) ' ihash theta thetas ', ihash, theta, 
     $           thetas(ihash)
            call errquit('xlm_n_m_l: duplicate theta hash',0)
         endif
      else
         call xlm_cache_theta(theta)
      endif
      itheta = ithetas(ihash)
c
      indp = 1
      indm = 1
      do l = 0, lmax
c
c     Add extra weights
c
         do m = -l, l
            qq(m,l) = q(m,l) * twobc(m,l)
         enddo
         qq(0,l) = qq(0,l) * 0.5d0
c
c     Apply Z rotation
c
         do m = 1, l
            qqp = cosmphi(m,iphi)*qq(m,l) + 
     $           phisign*sinmphi(m,iphi)*qq(-m,l)
            qqm = cosmphi(m,iphi)*qq(-m,l) -
     $           phisign*sinmphi(m,iphi)*qq(m,l)
            qq( m,l) = qqp
            qq(-m,l) = qqm
         enddo
c
c     Apply Y rotation
c
         do m = -l, l
            tmp(m) = qq(m,l)*rqn(m,l) ! Precondition
            qq(m,l) = 0.0d0
         enddo
         if (thetasign .gt. 0.0d0) then
            do mp = 0, l
               do m = 0, l
                  qq(m,l) = qq(m,l) + cdp(indp+m,itheta)*tmp(mp)
               enddo
               indp = indp + l + 1
            enddo
            do mp = -l,-1
               indm = indm + l
               do m = -l,-1
                  qq(m,l) = qq(m,l) + cdm(indm+m,itheta)*tmp(mp)
               enddo
            enddo
         else
            do m = 0, l
               do mp = 0, l
                  qq(m,l) = qq(m,l) + cdp(indp+mp,itheta)*tmp(mp)
               enddo
               indp = indp + l + 1
            enddo
            do m = -l,-1
               indm = indm + l
               do mp = -l,-1
                  qq(m,l) = qq(m,l) + cdm(indm+mp,itheta)*tmp(mp)
               enddo
            enddo
         endif
         do m = -l, l
            qq(m,l) = qq(m,l)*rtwobc(m,l)*qn(m,l) ! Remove precond & weights
         enddo
         qq(0,l) = qq(0,l) * 2.0d0
      enddo
c
c     Translate 0,0,r
c
      call xlm_multipole_to_local_z(lmax, thetasign*r, p, lp, qq, maxl)
c
c     Rotate the potentials back
c     
      indp = 1
      indm = 1
      do l = 0, lmax
c
c     Apply Y rotation
c
         do m = -l, l
            tmp(m) = p(m,l)*rqn(m,l) ! Precondition
            p(m,l) = 0.0d0
         enddo
         if (thetasign .lt. 0.0d0) then
            do mp = 0, l
               do m = 0, l
                  p(m,l) = p(m,l) + cdp(indp+m,itheta)*tmp(mp)
               enddo
               indp = indp + l + 1
            enddo
            do mp = -l,-1
               indm = indm + l
               do m = -l,-1
                  p(m,l) = p(m,l) + cdm(indm+m,itheta)*tmp(mp)
               enddo
            enddo
         else
            do m = 0, l
               do mp = 0, l
                  p(m,l) = p(m,l) + cdp(indp+mp,itheta)*tmp(mp)
               enddo
               indp = indp + l + 1
            enddo
            do m = -l,-1
               indm = indm + l
               do mp = -l,-1
                  p(m,l) = p(m,l) + cdm(indm+mp,itheta)*tmp(mp)
               enddo
            enddo
         endif
         do m = -l, l
            p(m,l) = p(m,l)*qn(m,l) ! Undo preconditioning
         enddo
c
c     Apply Z rotation
c
         do m = 1, l
            qqp = cosmphi(m,iphi)*p(m,l) - 
     $           phisign*sinmphi(m,iphi)*p(-m,l)
            qqm = cosmphi(m,iphi)*p(-m,l) +
     $           phisign*sinmphi(m,iphi)*p(m,l)
            p( m,l) = qqp
            p(-m,l) = qqm
         enddo
      enddo
c
      end
      subroutine xlm_cache_phi(phi)
      implicit none
#include "xlm.fh"
      double precision phi
c
c     Add to the list of phis with cached rotation info
c
      integer m, ihash
c
      if (numphi.lt.0 .or. numphi.ge.maxnumphi) call errquit
     $     ('xlm_cache_phi: numphi is invalid ', numphi)
      numphi = numphi + 1
c
      ihash = int(abs(phi*1000.0d0))
      if (iphis(ihash) .gt. 0) call errquit
     $     ('xlm_cache_phi: duplicate hash ', ihash)
      iphis(ihash) = numphi
      phis(ihash) = abs(phi)
c
      do m = 1, maxl
         cosmphi(m,numphi) = cos(dble(m)*abs(phi))
         sinmphi(m,numphi) = sin(dble(m)*abs(phi))
      end do
c
*      write(6,1) abs(phi), ihash, numphi
* 1    format(' xlm_cache_phi: added phi=',f9.6,
*     $     '   ihash=',i5,'   iphi=',i4)
c
      end
      subroutine xlm_cache_theta(theta)
      implicit none
#include "xlm.fh"
      double precision theta
c
c     Add to the list of thetas with cached rotation info
c
      integer l, ihash, indp, indm
c
      if (numtheta.lt.0 .or. numtheta.ge.maxnumtheta) call errquit
     $     ('xlm_cache_theta: numtheta is invalid ', numtheta)
      numtheta = numtheta + 1
c
      ihash = int(abs(theta*1000.0d0))
      if (ithetas(ihash) .gt. 0) call errquit
     $     ('xlm_cache_theta: duplicate hash ', ihash)
      ithetas(ihash) = numtheta
      thetas(ihash) = abs(theta)
c
      indp = 1
      indm = 1
      do l = 0, maxl
         call xlm_y_rotation_matrix(l, abs(theta),
     $        cdp(indp,numtheta), cdm(indm,numtheta))
         indp = indp + (l+1)*(l+1)
         indm = indm + l*l
      end do
      if (indm-1 .ne. dmdim) call errquit
     $     ('xlm_theta_cache: bad indm ', indm)
      if (indp-1 .ne. dpdim) call errquit
     $     ('xlm_theta_cache: bad indp ', indp)
c
*      write(6,1) abs(theta), ihash, numtheta
* 1    format(' xlm_cache_theta: added theta=',f9.6,
*     $     '   ihash=',i5,'   itheta=',i4)
c
      end

         
       
