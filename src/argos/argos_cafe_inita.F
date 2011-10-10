      subroutine argos_cafe_inita(nat,nap,nqt,nqp)
c
c     argos_cafe_init initializes atom type and charge type information
c
c     in:  nat, integer number of atom types
c          nap, integer number of parameters per atom type
c          nqt, integer number of charge types
c          nqp, integer number of parameters per charge type
c
c     this routine allocates the memory necessary to store all
c     atom type based force field parameters
c
      implicit none
c
      integer nat,nap,nqt,nqp
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      mat=max(1,nat)
      map=max(1,nap)
      mqt=max(1,nqt)
      mqp=max(1,nqp)
c
c     allocate memory
c
c     for the 4 van der Waals parameters (c6,c12,c6t,c12t) 1 : set 1
c                                                          2 : set 2
c                                                          3 : set 3
c                                                          4 : derivative
c                                                          5 : set 1
c                                                          6 : set 2
c
      if(.not.ma_push_get(mt_dbl,mat*mat*map*mset,'vdw',l_vdw,i_vdw))
     + call md_abort('Allocation failed for vdw',0)
c
c     for atomic masses 1 : mass set 1
c                       2 : mass set 2
c                       3 : mass set 3
c                       4 : mass derivative
c                       5 : mass set 1
c                       6 : mass set 2
c
      if(.not.ma_push_get(mt_dbl,mat*mset,'mas',l_mas,i_mas))
     + call md_abort('Allocation failed for mas',0)
c
c     for atom type names 1 : name set 1
c                         2 : name set 2
c                         3 : name set 3
c
      if(.not.ma_push_get(mt_byte,6*mat*nparms,'nam',l_nam,i_nam))
     + call md_abort('Allocation failed for nam',0)
c
c     for atomic numbers  1 : atomic number set 1
c                         2 : atomic number set 2
c                         3 : atomic number set 3
c
      if(.not.ma_push_get(mt_int,mat*nparms,'num',l_num,i_num))
     + call md_abort('Allocation failed for num',0)
c
c     for special type    1 : type in set 1
c                         2 : type in set 2
c                         3 : type in set 3
c
      if(.not.ma_push_get(mt_int,mat*nparms,'typ',l_typ,i_typ))
     + call md_abort('Allocation failed for typ',0)
c
      if(.not.ma_push_get(mt_int,mat*mat,'ias',l_ias,i_ias))
     + call md_abort('Allocation failed for ias',0)
c
c     for charge types
c
      if(.not.ma_push_get(mt_dbl,mqt*mqp*mset,'chg',l_chg,i_chg))
     + call md_abort('Allocation failed for chg',0)
c
      return
      end
c $Id$
