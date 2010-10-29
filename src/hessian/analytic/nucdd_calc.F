c
c $Id$
c
      subroutine nucdd_calc(c,zan,nat,dde)
c
c  This is a modified HONDO routine (hnd_nucdx) to calculate
c  the second derivative contribution of the nuclear term to
c  the hessian (dde).
c
      implicit none
c
#include "global.fh"
c
      integer nat                       ! [input] number of atoms
      double precision zan(nat)         ! [input] nuclear charges
      double precision c(3,nat)         ! [input] cartesian coordinates
      double precision dde(3,nat,3,nat) ! [input/output] hessian matrix
c
      integer iat, jat, ixyz, myid, nprocs
      double precision rij, dum, dum1, zero, three
c
      data zero   /0.0d+00/
      data three  /3.0d+00/
c
      if(nat.eq.1) go to 200
c
c     ----- calculate forces -----
c
      myid = ga_nodeid()  ! set up static parallel constructs
      nprocs = ga_nnodes()
c
      do 140 iat=myid+2,nat,nprocs   ! go parallel
      do 140 jat=1,iat-1
c
        rij=zero
        do 110 ixyz=1,3
          rij=rij+(c(ixyz,iat)-c(ixyz,jat))**2
  110   continue
        dum = zan(iat)*zan(jat)/(rij*sqrt(rij))
c
c     ----- ddedii , ddedij , ddedjj -----
c
        dum=dum/rij
        do 130 ixyz=1,3
          dum1=(three*(c(ixyz,iat)-c(ixyz,jat))**2-rij)*dum
          dde(ixyz,iat,ixyz,iat)=dde(ixyz,iat,ixyz,iat)+dum1
          dde(ixyz,iat,ixyz,jat)=dde(ixyz,iat,ixyz,jat)-dum1
          dde(ixyz,jat,ixyz,iat)=dde(ixyz,jat,ixyz,iat)-dum1
          dde(ixyz,jat,ixyz,jat)=dde(ixyz,jat,ixyz,jat)+dum1
  130   continue
        dum1=(three*(c(1,iat)-c(1,jat))*(c(2,iat)-c(2,jat)))*dum
        dde(1,iat,2,iat)=dde(1,iat,2,iat)+dum1
        dde(2,iat,1,iat)=dde(2,iat,1,iat)+dum1
        dde(1,iat,2,jat)=dde(1,iat,2,jat)-dum1
        dde(2,iat,1,jat)=dde(2,iat,1,jat)-dum1
        dde(1,jat,2,jat)=dde(1,jat,2,jat)+dum1
        dde(2,jat,1,jat)=dde(2,jat,1,jat)+dum1
c
        if (iat.ne.jat) then
          dde(1,jat,2,iat)=dde(1,jat,2,iat)-dum1
          dde(2,jat,1,iat)=dde(2,jat,1,iat)-dum1
        endif     
c
        dum1=(three*(c(1,iat)-c(1,jat))*(c(3,iat)-c(3,jat)))*dum
        dde(1,iat,3,iat)=dde(1,iat,3,iat)+dum1
        dde(3,iat,1,iat)=dde(3,iat,1,iat)+dum1
        dde(1,iat,3,jat)=dde(1,iat,3,jat)-dum1
        dde(3,iat,1,jat)=dde(3,iat,1,jat)-dum1
        dde(1,jat,3,jat)=dde(1,jat,3,jat)+dum1
        dde(3,jat,1,jat)=dde(3,jat,1,jat)+dum1
c
        if (iat.ne.jat) then
          dde(1,jat,3,iat)=dde(1,jat,3,iat)-dum1
          dde(3,jat,1,iat)=dde(3,jat,1,iat)-dum1
        endif     
c
        dum1=(three*(c(2,iat)-c(2,jat))*(c(3,iat)-c(3,jat)))*dum
        dde(2,iat,3,iat)=dde(2,iat,3,iat)+dum1
        dde(3,iat,2,iat)=dde(3,iat,2,iat)+dum1
        dde(2,iat,3,jat)=dde(2,iat,3,jat)-dum1
        dde(3,iat,2,jat)=dde(3,iat,2,jat)-dum1
        dde(2,jat,3,jat)=dde(2,jat,3,jat)+dum1
        dde(3,jat,2,jat)=dde(3,jat,2,jat)+dum1
c
        if (iat.ne.jat) then
          dde(2,jat,3,iat)=dde(2,jat,3,iat)-dum1
          dde(3,jat,2,iat)=dde(3,jat,2,iat)-dum1
        endif     
c
  140 continue
c
  200 continue
c
      return
      end
