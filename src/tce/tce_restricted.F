      subroutine tce_restricted_2(a1b,a2b,b1b,b2b)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer a1b,a2b
      integer b1b,b2b
c
c     Restricted
c
      if (restricted.and.
     1   (int_mb(k_spin+a1b-1)
     2   +int_mb(k_spin+a2b-1).eq.4)) then
        b1b = int_mb(k_alpha+a1b-1)
        b2b = int_mb(k_alpha+a2b-1)
      else
c
c     Unrestricted
c
       b1b = a1b
       b2b = a2b
      endif
      return
      end
      subroutine tce_restricted_4(a1b,a2b,a3b,a4b,
     1                            b1b,b2b,b3b,b4b)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer a1b,a2b,a3b,a4b
      integer b1b,b2b,b3b,b4b
c
c     Restricted
c
      if (restricted.and.
     1   (int_mb(k_spin+a1b-1)
     2   +int_mb(k_spin+a2b-1)
     3   +int_mb(k_spin+a3b-1)
     4   +int_mb(k_spin+a4b-1).eq.8)) then
        b1b = int_mb(k_alpha+a1b-1)
        b2b = int_mb(k_alpha+a2b-1)
        b3b = int_mb(k_alpha+a3b-1)
        b4b = int_mb(k_alpha+a4b-1)
      else
c
c     Unrestricted
c
       b1b = a1b
       b2b = a2b
       b3b = a3b
       b4b = a4b
      endif
      return
      end
      subroutine tce_restricted_6(a1b,a2b,a3b,a4b,a5b,a6b,
     1                            b1b,b2b,b3b,b4b,b5b,b6b)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer a1b,a2b,a3b,a4b,a5b,a6b
      integer b1b,b2b,b3b,b4b,b5b,b6b
c
c     Restricted
c
      if (restricted.and.
     1   (int_mb(k_spin+a1b-1)
     2   +int_mb(k_spin+a2b-1)
     3   +int_mb(k_spin+a3b-1)
     4   +int_mb(k_spin+a4b-1)
     5   +int_mb(k_spin+a5b-1)
     6   +int_mb(k_spin+a6b-1).eq.12)) then
        b1b = int_mb(k_alpha+a1b-1)
        b2b = int_mb(k_alpha+a2b-1)
        b3b = int_mb(k_alpha+a3b-1)
        b4b = int_mb(k_alpha+a4b-1)
        b5b = int_mb(k_alpha+a5b-1)
        b6b = int_mb(k_alpha+a6b-1)
      else
c
c     Unrestricted
c
       b1b = a1b
       b2b = a2b
       b3b = a3b
       b4b = a4b
       b5b = a5b
       b6b = a6b
      endif
      return
      end
      subroutine tce_restricted_8(a1b,a2b,a3b,a4b,a5b,a6b,a7b,a8b,
     1                            b1b,b2b,b3b,b4b,b5b,b6b,b7b,b8b)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer a1b,a2b,a3b,a4b,a5b,a6b,a7b,a8b
      integer b1b,b2b,b3b,b4b,b5b,b6b,b7b,b8b
c
c     Restricted
c
      if (restricted.and.
     1   (int_mb(k_spin+a1b-1)
     2   +int_mb(k_spin+a2b-1)
     3   +int_mb(k_spin+a3b-1)
     4   +int_mb(k_spin+a4b-1)
     5   +int_mb(k_spin+a5b-1)
     6   +int_mb(k_spin+a6b-1)
     7   +int_mb(k_spin+a7b-1)
     8   +int_mb(k_spin+a8b-1).eq.16)) then
        b1b = int_mb(k_alpha+a1b-1)
        b2b = int_mb(k_alpha+a2b-1)
        b3b = int_mb(k_alpha+a3b-1)
        b4b = int_mb(k_alpha+a4b-1)
        b5b = int_mb(k_alpha+a5b-1)
        b6b = int_mb(k_alpha+a6b-1)
        b7b = int_mb(k_alpha+a7b-1)
        b8b = int_mb(k_alpha+a8b-1)
      else
c
c     Unrestricted
c
       b1b = a1b
       b2b = a2b
       b3b = a3b
       b4b = a4b
       b5b = a5b
       b6b = a6b
       b7b = a7b
       b8b = a8b
      endif
      return
      end
