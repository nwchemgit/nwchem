      Subroutine mull_pop( geom, basis, iga_dens,iga_s, string)

C$Id: mull_pop.F,v 1.14 1998-10-14 23:43:12 windus Exp $
      Implicit none
      integer geom,basis
      integer iga_s              ! overlap   GA handle
      integer iga_dens           ! dens. mat GA handle
      integer iga_ps             ! product   GA handle
      character*(*) string      ! Identifier for ecce output
                                ! total/open/alpha/beta/mp2/ccsd/...

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
#include "cscfps.fh"
#include "mafdecls.fh"

      if (oscfps) call pstat_on(ps_mull)
      

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
     &     dbl_mb(iqatom),dbl_mb(iqshell),dbl_mb(iPSmat), string)
      
      call ga_sync

      if(.not.MA_Pop_Stack(lqshell))
     &     call errquit(' exiting in mull_pop',33)
      if(.not.MA_Pop_Stack(lqatom))
     &     call errquit(' exiting in mull_pop',34)
      if(.not.MA_Pop_Stack(lPSmat))
     &     call errquit(' exiting in mull_pop',35)
      status= ga_destroy(iga_PS)
c
      if (oscfps) call pstat_off(ps_mull)
c
      return
      end
c
c
c
      Subroutine mull_calc(basis, natoms, nshells,max_at_bf2,iga_PS,
     &     qatom,qshell,PSmat, string)

      Implicit none
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer basis
      integer natoms,nshells
      integer iga_PS             ! product   GA handle
      integer ifirst,ilast,nbf_at,max_at_bf2
      integer ish1,ish2,ish,nn,iat,mu
      integer me,nproc, geom
      double precision psmu, coord(3), qnuc
      double precision qatom(natoms),qshell(nshells),PSmat(max_at_bf2)
      character*2 symbol
      character*16 tag
      character*32 element
      character*80 buf
      character*(*) string
      integer lbuf
      integer atn
c     
      me=ga_nodeid()
      nproc=ga_nnodes()
      
      call dfill(natoms,0.D0,qatom,1)
      call dfill(nshells,0.D0,qshell,1)

      if (.not. bas_geom(basis, geom)) call errquit
     $     ('mull_pop: bas_geom failed',basis)

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
    1   format(/'    Atom       Charge   Shell Charges')
        write(luout,11)
 11     format( ' -----------   ------   -----------------------------',
     $       '--------------------------')
        do iat=1,natoms
          if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &         call errquit(' exiting in mull_pop',4)
c
          if (.not. geom_cent_get(geom, iat, tag, coord, qnuc)) 
     $         call errquit('mull_pop: geom_cent_tag failed',0)
c
          if (.not. geom_tag_to_element(tag, symbol, element, atn))
     $         symbol = 'X'
          if (ish2.ge.ish1) then
             write(LuOut,2) iat,symbol,nint(qnuc),qatom(iat),
     &            (qshell(ish),ish=ish1,ish2)
 2           format(1x,i4,1x,a2,i4,1x,f8.2,1x,25(f6.2))
          endif
        enddo
      endif
c
      buf = ' '
      write(buf,'(a,a)') string(1:inp_strlen(string)), 
     $     ' mulliken atomic charges'
      lbuf = inp_strlen(buf)
c
c     Calculate total atomic charges for ecce
c
      do iat=1,natoms
        if (.not. geom_cent_get(geom, iat, tag, coord, qnuc))
     $       call errquit('mull_pop: geom_cent_tag failed',0)
        qatom(iat) = qnuc - qatom(iat)
      enddo
c
      call ecce_print1(buf(1:lbuf), mt_dbl, qatom, natoms)
      buf = ' '
      write(buf,'(a,a)') string(1:inp_strlen(string)), 
     $     ' mulliken shell charges'
      lbuf = inp_strlen(buf)
      call ecce_print1(buf(1:lbuf), mt_dbl, qshell,nshells)
      call ga_sync

      return
      end
