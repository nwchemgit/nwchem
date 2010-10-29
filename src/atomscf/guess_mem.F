      subroutine guess_mem(memscr)
C$Id$
      implicit none
#include "cguess.fh"
c
c (atw) Modified to calculate max. memory for atom blocks
c       only. 7/3/94
c
      integer memscr
c
c..   calculate memory requirements for atomic guess routines
c..
      integer ntr, nsq
      integer i10, ipcap, iqcap, ifc, ifo, is, iu, it
      integer ih, idc, idos, idt, idold, iss, ic, icopn, ismin
      integer iqmin, itransf, icc
c..
c..    core partitioning
c..

c
c
      ntr = nb*(nb+1)/2
      nsq = nb * nb
c
c
      i10 = 1
      ipcap = i10 + ntr
      iqcap = ipcap + ntr
      ifc = iqcap + ntr
      ifo = ifc + ntr
      is = ifo + ntr
      iu = is + ntr
      it = iu + ntr
      ih = it + ntr
      idc = ih + ntr
      idos = idc + ntr
      idt = idos + ntr
      idold = idt + ntr
      iss = idold + ntr
      ic = iss + ntr
c
      icopn = ic + nsq
      ismin = icopn + nsq
      iqmin = ismin + nb * no
      itransf = iqmin + nb * no
      icc = itransf + nsq
      memscr = icc + nsq - 1
c
c     NOTE: later inserts required for pseudopotentials
c
      end
