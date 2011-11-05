      LOGICAL FUNCTION IS_ACTIVE_1(g1b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER g1b
      is_active_1 = .false.
      if(int_mb(k_active+g1b-1).eq.1) is_active_1 = .true.
      RETURN
      END
      LOGICAL FUNCTION IS_ACTIVE_2(g1b,g2b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER g1b
      INTEGER g2b
      is_active_2 = .false.
      if((int_mb(k_active+g1b-1).eq.1).and.
     &   (int_mb(k_active+g2b-1).eq.1)) is_active_2 = .true.
      RETURN
      END
      LOGICAL FUNCTION IS_ACTIVE_3(g1b,g2b,g3b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER g1b
      INTEGER g2b
      INTEGER g3b
      is_active_3 = .false.
      if((int_mb(k_active+g1b-1).eq.1).and.
     &   (int_mb(k_active+g2b-1).eq.1).and.
     &   (int_mb(k_active+g3b-1).eq.1)) is_active_3 = .true. 
      RETURN
      END
      LOGICAL FUNCTION IS_ACTIVE_4(g1b,g2b,g3b,g4b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER g1b
      INTEGER g2b
      INTEGER g3b
      INTEGER g4b
      is_active_4 = .false.
      if((int_mb(k_active+g1b-1).eq.1).and.
     &   (int_mb(k_active+g2b-1).eq.1).and.
     &   (int_mb(k_active+g3b-1).eq.1).and.
     &   (int_mb(k_active+g4b-1).eq.1)) is_active_4 = .true.
      RETURN
      END
      LOGICAL FUNCTION IS_ACTIVE_4_O(g1b,g2b,g3b,g4b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER g1b
      INTEGER g2b
      INTEGER g3b
      INTEGER g4b
      is_active_4_o = .false.
      if((int_mb(k_active_o+g1b-1).eq.1).and.
     &   (int_mb(k_active_o+g2b-1).eq.1).and.
     &   (int_mb(k_active_o+g3b-1).eq.1).and.
     &   (int_mb(k_active_o+g4b-1).eq.1)) is_active_4_o = .true.
      RETURN
      END
      LOGICAL FUNCTION ONE_OF_TWO_ACT(g1b,g2b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER g1b
      INTEGER g2b
      one_of_two_act = .false.
      if((int_mb(k_active+g1b-1).eq.1).or.
     &   (int_mb(k_active+g2b-1).eq.1)) one_of_two_act = .true.
      RETURN
      END
      LOGICAL FUNCTION is_active_t(h1,h2,h3,p1,p2,p3)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1,h2,h3
      INTEGER p1,p2,p3
      INTEGER ahole,aparticle
      is_active_t = .false.
      ahole=int_mb(k_active+h1-1)+int_mb(k_active+h2-1)
     &     +int_mb(k_active+h3-1)
      aparticle=int_mb(k_active+p1-1)+int_mb(k_active+p2-1)
     &      +int_mb(k_active+p3-1)
      if((ahole.ge.2).and.(aparticle.ge.2)) is_active_t = .true.
      RETURN
      END
      LOGICAL FUNCTION is_active_t_iii(h1,h2,h3,p1,p2,p3)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1,h2,h3
      INTEGER p1,p2,p3
      INTEGER ahole,aparticle
      is_active_t_iii = .false.
      ahole=int_mb(k_active+h1-1)+int_mb(k_active+h2-1)
     &     +int_mb(k_active+h3-1)
      aparticle=int_mb(k_active+p1-1)+int_mb(k_active+p2-1)
     &      +int_mb(k_active+p3-1)
      if((ahole.eq.3).and.(aparticle.eq.3)) is_active_t_iii = .true.
      RETURN
      END

c $Id$
