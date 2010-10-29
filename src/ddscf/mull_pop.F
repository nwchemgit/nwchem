      Subroutine mull_pop( geom, basis, iga_dens,iga_s, string)

C$Id$
      Implicit none
#include "errquit.fh"
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
     &     call errquit(' exiting in mull_pop',0, GEOM_ERR)
      if( .not. bas_numcont(basis,nshells) )
     &     call errquit(' exiting in mull_pop',1, BASIS_ERR)
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in mull_op',1, BASIS_ERR)
      max_at_bf2 = 0
      do iatom = 1, natoms
        if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $       call errquit('mul_pop: bas_ce2bfr failed', iatom,
     &       BASIS_ERR)
        max_at_bf2 = max(max_at_bf2, ihi-ilo+1)
      enddo
      max_at_bf2 = max_at_bf2*max_at_bf2

      if(.not.MA_Push_Get(mt_dbl,max_at_bf2,'PS',lPSmat,iPSmat))
     &     call errquit(' exiting in mull_pop: insuff stack',21, MA_ERR)
      if(.not.MA_Push_Get(mt_dbl,natoms,'q atom',lqatom,iqatom))
     &     call errquit(' exiting in mull_pop: insuff stack',22, MA_ERR)
      if(.not.MA_Push_Get(mt_dbl,nshells,'q shell',lqshell,iqshell))
     &     call errquit(' exiting in mull_pop: insuff stack',3, MA_ERR)

      iga_PS=ga_create_atom_blocked(geom, basis, 'PS product')

      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,
     &     iga_dens,iga_s,0.d0,iga_PS)
      call mull_calc(basis,natoms, nshells,max_at_bf2,iga_PS,
     &     dbl_mb(iqatom),dbl_mb(iqshell),dbl_mb(iPSmat), string)
      
      call ga_sync

      if(.not.MA_Pop_Stack(lqshell))
     &     call errquit(' exiting in mull_pop',33, MA_ERR)
      if(.not.MA_Pop_Stack(lqatom))
     &     call errquit(' exiting in mull_pop',34, MA_ERR)
      if(.not.MA_Pop_Stack(lPSmat))
     &     call errquit(' exiting in mull_pop',35, MA_ERR)
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
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      integer basis
      integer natoms,nshells
      integer iga_PS             ! product   GA handle
      integer ifirst,ilast,nbf_at,max_at_bf2
      integer ish1,ish2,ish,nn,iat,mu
      integer me,nproc, geom
      double precision psmu, coord(3), qnuc, pop
      double precision qatom(natoms),qshell(nshells),PSmat(max_at_bf2)
      character*2 symbol
      character*16 tag
      character*32 element
      character*80 buf
      character*(*) string
      integer lbuf
      integer atn
      integer nbf, k_tags, l_tags, g_nbf, i, m, offset, nbf_acc, g_pop
      logical oprint
c     
      me=ga_nodeid()
      nproc=ga_nnodes()
      oprint = util_print('Mulliken AO', print_high)
      
      call dfill(natoms,0.D0,qatom,1)
      call dfill(nshells,0.D0,qshell,1)

      if (.not. bas_geom(basis, geom)) call errquit
     $     ('mull_pop: bas_geom failed',basis, BASIS_ERR)

       if (.not.bas_numbf(basis,nbf))
     $     call errquit('bas_numbf failed',0, BASIS_ERR)

      if(.not.ga_create(mt_int,nbf,1,'atom nbf',0,0,g_nbf))
     $        call errquit('mull_pop ga_create failed', 0, GA_ERR)
      call ga_zero(g_nbf)
 
      if(.not.ga_create(mt_dbl,nbf,natoms,'ao pop',0,0,g_pop))
     $        call errquit('mull_pop ga_create failed', 0, GA_ERR)
      call ga_zero(g_pop)
 
      if (.not. ma_push_get(mt_byte,nbf*16,'bf tags',l_tags,k_tags))
     $    call errquit(' exiting in mull_pop: insuff stack',3, MA_ERR)

      call bas_vec_info(basis, byte_mb(k_tags))

 3    format(1x,(i5,2x,f12.6,4x,16a1,4x))
 4    format(3x,a4,4x,a10,2x,a13)
      if (me.eq.0 .and. oprint) then
        write(6,*)
        write(6,4)'Bfn.','Population','Atom+Function'
        write(6,*)'  ---------------------------------'
      endif

      call ga_sync()

      do iat=me+1,natoms,nproc
        if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &       call errquit(' exiting in mull_pop',4, BASIS_ERR)
        call get_atom_block(iga_PS, basis,
     $       iat, iat, PSmat, nbf_at, nbf_at)
        call ga_put(g_nbf,iat,iat,1,1,nbf_at,1)

        mu=0
        do ish=ish1,ish2
          if (.not. bas_cn2bfr(basis,ish,ifirst,ilast))
     &         call errquit(' exiting in mull_pop.',5, BASIS_ERR)
          do nn=ifirst,ilast
            mu=mu+1
            psmu=PSmat((mu-1)*nbf_at+mu)
            
            call ga_put(g_pop,mu,mu,iat,iat,psmu,1)
 
            qshell(ish)=qshell(ish)+psmu
          enddo
          qatom(iat)=qatom(iat)+qshell(ish)
        enddo
      enddo

      call ga_sync()
     
      if (me.eq.0 .and. oprint) then
      offset = 0
      do iat = 1,natoms
        call ga_get (g_nbf,iat,iat,1,1,nbf_at,1)
        do i = 1,nbf_at
          call ga_get(g_pop,i,i,iat,iat,pop,1)
          write(6,3) i+offset,pop,
     &     (byte_mb(k_tags+(i-1+offset)*16+m),m=0,15)
        enddo
          offset = offset + nbf_at
       enddo
       endif

      if(.not.ga_destroy(g_pop))
     $        call errquit('mull_pop ga_destroy failed', 0, GA_ERR)

      if(.not.ga_destroy(g_nbf))
     $        call errquit('mull_pop ga_destroy failed', 0, GA_ERR)

      if(.not.MA_Pop_Stack(l_tags))
     &     call errquit(' exiting in mull_pop',32, MA_ERR)

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
     &         call errquit(' exiting in mull_pop',4, BASIS_ERR)
c
          if (.not. geom_cent_get(geom, iat, tag, coord, qnuc)) 
     $         call errquit('mull_pop: geom_cent_tag failed',0,
     &       GEOM_ERR)
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
     $       call errquit('mull_pop: geom_cent_tag failed',0, GEOM_ERR)
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
