*
* $Id: stpr_gen_aux.F,v 1.1 1995-04-24 08:04:24 d3e129 Exp $
* 
* this file contains auxilary routines to stpr_gen functions
*
* current routines: 
* stpr_fd_upd_hess     ! computes finite difference either central or forward
* stpr_wrt_fd_from_sq  ! writes hessian to file
c
      subroutine stpr_fd_upd_hess(hess,gradm,gradp,s_delta,delta,nat,
     &    iatom_t,ixyz_t)
      implicit none
#include "stdio.fh"
c::passed
      integer nat
      integer iatom_t, ixyz_t
      double precision hess(3,nat,3,nat)
      double precision gradm(3,nat)
      double precision gradp(3,nat)
      double precision delta, s_delta
c::local
      integer iatom,ixyz
      double precision rdelta, value
c
c 
c finite difference  [g(x+delta) - g(x-delta)]/(s_delta*delta) (s_delta = 2.0)
c central difference [g(x+delta) - g(x)]/(s_delta*delta) (s_delta = 1.0)
c 
c
      rdelta = 1.0d00/(s_delta*delta)
      do 00100 iatom = 1,nat
        do 00200 ixyz = 1,3
          value = rdelta*(gradp(ixyz,iatom)-gradm(ixyz,iatom))
          hess(ixyz_t,iatom_t,ixyz,iatom) = value
          hess(ixyz,iatom,ixyz_t,iatom_t) = value
00200   continue
00100 continue
c
      end
      subroutine stpr_wrt_fd_from_sq(hess,rank_hess,filename)
      implicit none
c
#include "stdio.fh"
c
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
      character*(*) filename
c
      logical does_it_exist
c
      integer i, j, lu
c
      lu = 66
      does_it_exist = .false.
      inquire(file=filename,exist=does_it_exist)
      if (does_it_exist)
     &    write(luout,*)
     &    ' stpr_wrt_fd_from_sq: overwrite of existing file:',
     &    filename
      open(unit=lu,file=filename,
     &    form='formatted',
     &    access='sequential',
     &    status='unknown')
c
      do 00100 i = 1,rank_hess
        do 00200 j = 1,i
          write(lu,10000)hess(i,j)
00200   continue
00100 continue
c
10000 format(1x,1pd20.10)
c
      close(unit=lu,status='keep')
c
      end
