      subroutine int_acc_std()
C $Id: int_acc.F,v 1.5 1995-11-29 02:10:36 d3e129 Exp $
      implicit none
#include "apiP.fh"
      double precision val_def
      parameter(val_def = 1.0d-10)
c
      val_int_acc = val_def
      call int_sp_acc(val_int_acc)
      end
c......................................................................
      subroutine int_acc_high()
      implicit none
#include "apiP.fh"
      val_int_acc = 1.0d-20
      call int_sp_acc(val_int_acc)
      call texas_set_accy(val_int_acc,val_int_acc)
      end
c......................................................................
      subroutine int_acc_get(retval)
      implicit none
#include "apiP.fh"
      double precision retval
      retval = val_int_acc
      end
c......................................................................
      subroutine int_acc_set(setval)
      implicit none
#include "apiP.fh"
      double precision setval
      val_int_acc = setval
      call int_sp_acc(val_int_acc)
      call texas_set_accy(val_int_acc,val_int_acc)
      end
