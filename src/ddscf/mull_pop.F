      Subroutine mull_pop( geom, basis, iga_dens,iga_s)

C$Id: mull_pop.F,v 1.2 1996-05-15 17:41:41 og845 Exp $
      Implicit none
      integer geom,basis
      integer iga_s              ! overlap   GA handle
      integer iga_dens           ! dens. mat GA handle
      integer iga_ps             ! product   GA handle

      integer natoms,nshells
      integer lPSmat,iPSmat,lqatom,iqatom,lqshell,iqshell
      integer iatom,ilo,ihi,nbf,max_at_bf2
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical status

#include "bas.fh"
#include "geom.fh"
#include "global.fh"

#include "mafdecls.fh"


      

c*****************************************************************************

c     
      if(.not.geom_ncent(geom, natoms))
     &     call errquit(' exiting in mull_pop',0)
      if( .not. bas_numcont(basis,nshells) )
     &     call errquit(' exiting in mull_pop',1)
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in mull_op',1)
      max_at_bf2 = 0
      do iatom = 1, natoms
        if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $       call errquit('mul_pop: bas_ce2bfr failed', iatom)
        max_at_bf2 = max(max_at_bf2, ihi-ilo+1)
      enddo
      max_at_bf2 = max_at_bf2*max_at_bf2

      if(.not.MA_Push_Get(mt_dbl,max_at_bf2,'PS',lPSmat,iPSmat))
     &     call errquit(' exiting in mull_pop: insuff stack',21)
      if(.not.MA_Push_Get(mt_dbl,natoms,'q atom',lqatom,iqatom))
     &     call errquit(' exiting in mull_pop: insuff stack',22)
      if(.not.MA_Push_Get(mt_dbl,nshells,'q shell',lqshell,iqshell))
     &     call errquit(' exiting in mull_pop: insuff stack',3)

      iga_PS=ga_create_atom_blocked(geom, basis, 'PS product')

      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,
     &     iga_dens,iga_s,0.d0,iga_PS)
      call mull_calc(basis,natoms, nshells,max_at_bf2,iga_PS,
     &     dbl_mb(iqatom),dbl_mb(iqshell),dbl_mb(iPSmat))
      
      call ga_sync

      if(.not.MA_Pop_Stack(lqshell))
     &     call errquit(' exiting in mull_pop',33)
      if(.not.MA_Pop_Stack(lqatom))
     &     call errquit(' exiting in mull_pop',34)
      if(.not.MA_Pop_Stack(lPSmat))
     &     call errquit(' exiting in mull_pop',35)
      status= ga_destroy(iga_PS)
      return
      end
c
c
c
      Subroutine mull_calc(basis, natoms, nshells,max_at_bf2,iga_PS,
     &     qatom,qshell,PSmat)

      Implicit none
      integer basis
      integer natoms,nshells
      integer iga_PS             ! product   GA handle

      integer ifirst,ilast,nbf_at,max_at_bf2
      integer ish1,ish2,ish,nn,iat,mu
      integer me,nproc
      double precision psmu
      double precision qatom(natoms),qshell(nshells),PSmat(max_at_bf2)

#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "msgids.fh"

c     
      me=ga_nodeid()
      nproc=ga_nnodes()
      
      call dfill(natoms,0.D0,qatom,1)
      call dfill(nshells,0.D0,qshell,1)

      do iat=me+1,natoms,nproc
        if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &       call errquit(' exiting in mull_pop',4)
        call get_atom_block(iga_PS, basis,
     $       iat, iat, PSmat, nbf_at, nbf_at)
        mu=0
        do ish=ish1,ish2
          if (.not. bas_cn2bfr(basis,ish,ifirst,ilast))
     &         call errquit(' exiting in mull_pop.',5)
          do nn=ifirst,ilast
            mu=mu+1
            psmu=PSmat((mu-1)*nbf_at+mu)
            qshell(ish)=qshell(ish)+psmu
          enddo
          qatom(iat)=qatom(iat)+qshell(ish)
        enddo
      enddo
      call ga_sync
      call ga_dgop(Msg_Mull1,qatom,natoms,'+')
      call ga_dgop(Msg_Mull2,qshell,nshells,'+')
      if(me.eq.0) then
        write(LuOut,1)
    1   format(/' atom     Charge     Shell Charges')
        do iat=1,natoms
          if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &         call errquit(' exiting in mull_pop',4)
          write(LuOut,2) iat,qatom(iat),
     &         (qshell(ish),ish=ish1,ish2)
    2     format(i4,2x,f10.4,20(f10.4))
        enddo
      endif
      call ga_sync

      return
      end
