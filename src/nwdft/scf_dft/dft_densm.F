c
      subroutine dft_densm(g_dens,g_vecs,nbf,nmo,ipol,
     ,     geom,AO_bas_han,noc,ntotel,
     ,	evals,focc,ssmear,tdots,
     ,     iter,search0,
     .     fon, nel_fon,nmo_fon,ncore_fon,
     .                     spinset, rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "dftps.fh"
#include "util.fh"
#include "rtdb.fh" 
c
      integer geom
      integer AO_bas_han
      integer ipol                  ! no. of spin states
      integer g_dens(2)
      integer g_vecs(2)
      integer nbf,nmo(2)
      logical spinset
      integer noc(2)             ! no. of occupied levels
      integer ntotel,iter
      double precision focc(nbf*ipol)! occupation no.
      double precision evals(ipol*nbf)
      double precision ssmear       ! smearing sigma
      double precision tdots       ! total energy
      logical search0,fon
      double precision nel_fon(4)
      integer ncore_fon(2), nmo_fon(4)
      integer ilo,ihi,jlo,jhi,numi,numj,icount
c
      integer i,isp,g_tmp
      integer me, nproc
      double precision zero,toll,efermi(2)
      double precision cksum,ntot(2)
      integer ltmpm,itmpm,ipoint
      double precision sqrtpi,kbau,eps,rhfuhf
      double precision one,x,xb
      double precision util_erfc
      external util_erfc
      integer nocsave(2)
c
      logical dft_fon
      external dft_fon
c
      double precision undovl
      parameter(undovl=-20d0*2.3025d0)
      parameter (zero=0.d0,toll=1.d-9,one=1.d0,
     ,   kbau=1.d0,eps=1.d-4)
c
      integer rtdb         ! FA-02-22-11, for occupations
      integer switch_focc  ! FA-02-22-11, for occupations
      logical status       ! FA-02-22-11, for occupations
      double precision anoc(2),new_ntotel
      logical debug_fon
      logical no_avg_fon
      integer nocinit
      save nocinit,anoc
      data nocinit /1/
c
c     preliminaries
      sqrtpi=sqrt(acos(-1d0))
      me=ga_nodeid()
      nproc=ga_nnodes()
      rhfuhf=2d0
      nocsave(1)=-1
      nocsave(2)=-2
      if(ipol.eq.2) rhfuhf=1d0
c       
c     check degener
c       
        if (rtdb_get(rtdb, 'dft:debugfon', mt_log, 1,
     &     debug_fon)) then
          if (debug_fon .and. me.eq.0)
     &       write(luout,*) "fractional occupation (fon) is on ",fon
        end if
c
c       average fractional occupations (leading orbitals)
        if (fon) then
          status = dft_fon(rtdb,nmo_fon,ncore_fon,nel_fon,
     .       nbf,ntotel,focc,noc,ipol,me)
          do isp=1,ipol
            call dft_focdm(focc(1+(isp-1)*nbf),noc(isp),geom,
     &         AO_bas_han,nbf,g_vecs(isp),g_dens(isp),toll)
          enddo
          if(ssmear.lt.toll)return
        endif ! fon
c
      if(ssmear.lt.toll.or.iter.lt.-1) then
c
        if (odftps) call pstat_on(ps_dgemm)
c
c       fractional occupation by orbital
        switch_focc=0
        status=rtdb_get(rtdb,'focc:occup_switch', mt_int,1,switch_focc)
        if (switch_focc.eq.1 .and. status) then ! using specified occupations
          do isp=1,ipol
            call dft_frac_mo(rtdb,g_vecs(isp),focc,nbf,ipol,ntotel)
            call dft_focdm(focc(1+(isp-1)*nbf),noc(isp),geom,
     &         AO_bas_han,nbf,g_vecs(isp),g_dens(isp),toll)
          enddo  ! isp
          return
        else                    ! default occupations
          do isp=1,ipol
            call ga_dgemm('N', 'T', nbf, nbf, noc(isp),
     $         2d0/dble(ipol), g_vecs(isp),g_vecs(isp),
     $         zero, g_dens(isp))
          enddo
        endif                   ! switch_focc

        if (odftps) call pstat_off(ps_dgemm)
c
      else	  ! with smearing
c
        if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &     call errquit('dftdensm: failed to alloc tmpm',nbf, MA_ERR)
c
c     fon
c
        if(nocinit.eq.1) then
           if(fon) then
              anoc(1)=ncore_fon(1)+nel_fon(1)
              anoc(2)=0d0
              if(ipol.eq.2)anoc(2)=ncore_fon(2)+nel_fon(2)
              nocinit=0
           else
              anoc(1)=noc(1)
              anoc(2)=0d0
              if(ipol.eq.2) anoc(2)=noc(2)
              nocinit=1
           endif
        endif
        new_ntotel=rhfuhf*(anoc(1)+anoc(2))
c       
c       initialize ef
c       
        efermi(1)=0.d0
        efermi(2)=0.d0
        if (.not. ga_duplicate(g_dens, g_tmp, 'ga_temp'))
     &     call errquit('dft_densm: error creating ga',0,GA_ERR)
        call ga_zero(g_tmp)
        if(search0)then
          call dfill(nbf*ipol, 0.d0, focc, 1)
          if(spinset.and.ipol.eq.2) then
            nmo(2)=0
            call dft_zero(2,nbf,nmo,anoc(1),efermi(1),evals,
     ,         ssmear,toll,.true.)
            call dft_zero(2,nbf,nmo,anoc(2),efermi(2),evals(nbf+1),
     .         ssmear,toll,.true.)
            nmo(2)=nmo(1)
          else
cold            call dft_zero(ipol,nbf,nmo,ntotel,efermi(1),evals,ssmear,
            call dft_zero(ipol,nbf,nmo,new_ntotel,
     E            efermi(1),evals,ssmear,
     .         toll,spinset)
            efermi(2)=efermi(1)
          endif
        else
          if(spinset) then
            nocsave(1)=noc(1)
            if(ipol.eq.2) nocsave(2)=noc(2)
          endif
        endif
        do isp=1,ipol
          ntot(isp)=0.d0
          ipoint=(isp-1)*nbf
          do i = 1, nmo(isp)
            x=(evals(ipoint+i)-efermi(isp))/ssmear
            if(search0)
     *         focc(i+ipoint)= util_erfc(x)*0.5d0*rhfuhf
            if(focc(i+ipoint).gt.toll) then
              noc(isp)=i
            else
              focc(i+ipoint)=0.d0
            endif
          enddo
          call ga_zero(g_tmp)
          call ga_distribution(g_vecs(isp),
     .       ga_nodeid(), ilo, ihi, jlo, jhi)
          if (ilo.gt.0 .and. ilo.le.ihi) then
            numi =  ihi-ilo+1
            numj =  jhi-jlo+1
            if (numi.gt.0 .and. numj.gt.0) then
              if (.not. MA_Push_Get(MT_Dbl, numi*numj, 'matS',
     T           ltmpm, itmpm)) call errquit(
     T           'dftdensm: failed to alloc tmpm',0,0)
              call ga_get(g_vecs(isp),ilo,ihi,jlo,jhi,
     .           dbl_mb(itmpm),numi)
              icount=itmpm
              do i=jlo,jhi
                if(focc(i+ipoint).gt.toll) then
                  call dscal(numi, focc(i+ipoint), 
     D               DBL_MB(icount), 1)
                else
                  call dcopy(numi, 0d0,0, DBL_MB(icount), 1)
                endif
                icount=icount+numi
              enddo
              call ga_put(g_tmp,ilo,ihi,jlo,jhi,
     .           dbl_mb(itmpm),numi)
              if (.not.ma_pop_stack(ltmpm)) call errquit
     T           ('dft_densm: cannot pop stack',0, MA_ERR)
            endif
          endif 
          do i=1,nmo(isp)
            ntot(isp)=ntot(isp)+focc(i+ipoint)
          enddo
          call ga_sync()
          call ga_dgemm('N', 'T', nbf, nbf, noc(isp),1.d0,g_tmp,
     $       g_vecs(isp), zero, g_dens(isp))
        enddo                   ! end of loop over spins
        if (.not. ga_destroy(g_tmp)) call errquit
     &     ('dft_densm: could not destroy g_tmp', 0, GA_ERR)
c       
c       compute entropy correction to total e
c     
        if(ssmear.gt.1d-8) then
          cksum=0.d0
         do i=1,nmo(1)
           x=(evals(i)-efermi(1))/ssmear
           if((-x*x).gt.undovl)
     a        cksum=cksum+exp(-x*x)/sqrtpi
           if(ipol.eq.2) then
             xb=(evals(i+nbf)-efermi(2))/ssmear
             if((-xb*xb).gt.undovl)
     b          cksum=cksum+exp(-xb*xb)/sqrtpi
           endif
         enddo
         tdots=tdots-cksum*ssmear*0.5d0*rhfuhf
         if(ga_nodeid().eq.0 .and.
     $      util_print('general information',print_default) ) then 
           write(LuOut,1688) ssmear,
     ,        efermi(1),tdots,(i,ntot(i),i=1,ipol)
 1688      format(' sigma= ',d9.2,
     ,        ' eFermi= ',f8.5,' T*S= ',1p,e9.2,0p,
     ,        2(' nel(',i1,')=',f8.2)/) 
C          write(LuOut,1723) (focc(i),i=1,nbf*ipol)
C          1723       format(' Foccs ',(
         endif
       endif
c
       if(spinset) then
         if(search0) then
           noc(1)=nint(ntot(1)*ipol*0.5d0)
           if(ipol.eq.2) noc(2)=nint(ntot(2))
         else
           noc(1)=nocsave(1)
           noc(2)=nocsave(2)
         endif
       endif
c
      endif  ! smear
c
      return
      end
c
      subroutine dft_zero(ipol,nbf,nmo,
     ,     ntotel,efermi,evals,ssmear,toll,spinset)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer ipol
      integer nbf,nmo(2)
      double precision evals(nbf*ipol) 
      double precision ssmear
      double precision toll
cold      integer ntotel
      double precision ntotel
      double precision efermi
      logical spinset
c
      integer attempts,mx_attem
      double precision emax,emin,fmax,fmin,fmed,step
      parameter(mx_attem=200)
      external fzero
      double precision fzero

      attempts=0
C
C     determine extrema of search
C
      if (ipol.eq.1) then
c
c     closed shell
c
         efermi=evals(int(ntotel/2))
      else
c
c     open shell
c
         if(spinset) then
            efermi=evals(int(ntotel/2))
         else
            efermi=0.5d0*(evals(int(ntotel/2))+
     +           evals(nbf+int(ntotel/2)))
         endif
      endif
      step=max(2.d0*ssmear,1.d0)
 11   emax=efermi+step
      emin=efermi-step
      fmax=fzero(
     *     ipol,emax,evals,nbf,nmo,ntotel,ssmear)
      fmin=fzero(
     *     ipol,emin,evals,nbf,nmo,ntotel,ssmear)
      fmed=fzero(
     *     ipol,efermi,evals,nbf,nmo,ntotel,ssmear)
C      if(abs(fmed).lt.toll) goto
      if(fmax*fmin.gt.0.d0) then
         attempts=attempts+1
      if(attempts.gt.mx_attem) then
         if(ga_nodeid().eq.0) write(0,*) ' efermi emax emin ',
     ,        efermi,emax,emin
         if(ga_nodeid().eq.0) write(0,*) ' fmed fmax fmin ',
     ,        fmed,fmax,fmin
       call errquit('dft_zero: extrema:too many attempts',mx_attem,
     &       CALC_ERR)
      endif
         step=step+step
         goto 11
      endif
      attempts=0
c
c     bisection brute-force search
c
 1    continue
      attempts=attempts+1
      if(attempts.gt.mx_attem) then
#ifdef DEBUG
         if(ga_nodeid().eq.0) then
         write(LuOut,*) ' efermi emax emin ',
     ,        efermi,emax,emin
         write(LuOut,*) ' fmed fmax fmin ',fmed,fmax,fmin
         write(LuOut,*) 'ZERO could no reach required accuracy'
         endif
#endif
C         call errquit('dft_zero: too many attempts',mx_attem)     
         goto 1999
       endif
      if(fmax*fmed.gt.0) then
c	zero between emin and efermi
         emax=efermi
         fmax=fmed
      else	
c     zero between fmed and fmax
         emin=efermi
         fmin=fmed
      endif
C         if(me.eq.0) then
C            write(0,*) ' fmax = ',fmax,' fmin = ',fmin
C            write(0,*) ' emax = ',emax,' emin = ',emin
C            write(0,*) ' attempts ',attempts
C         endif
      if(abs(fmax-fmin).le.0.d0) then
         if(ga_nodeid().eq.0) then
            write(0,*) ' WARNING: fmax=fmin'
         endif
      endif
      if(attempts.lt.15.or.abs(fmax-fmin).lt.0.d0) then
      efermi=0.5d0*(emin+emax)
      else
         efermi=-fmin*(emax-emin)/(fmax-fmin)+emin
      endif
      
      fmed=fzero(
     *     ipol,efermi,evals,nbf,nmo,ntotel,ssmear)
      if(abs(fmed).gt.toll.and.abs(emin-emax).gt.toll*0.1d0) goto 1
1999  continue
#ifdef DEBUG
      if(ga_nodeid().eq.0) write(LuOut,*) ' fmed emin-emax ',abs(fmed),
     *abs(emin-emax)
#endif

      return
      end
c
      double precision function fzero(
     *     ipol,efermi,evals,nbf,nmo,ntotel,ssmear)
      implicit none
      integer ipol
      integer nbf,nmo(2)
      double precision efermi
      double precision evals(nbf*ipol)
cold      integer ntotel 
      double precision ntotel
      double precision ssmear
#include "msgids.fh"
c
      integer i
      double precision util_erfc
      external util_erfc
      double precision one,x,xb,rhfuhf,
     , fzeroa ,fzerob ,na,nb
      rhfuhf=2d0/ipol
      fzeroa=0d0
      fzerob=0d0
      nb=0d0
      do i=1,nmo(1)
         x=(evals(i)-efermi)/ssmear
         na=util_erfc(x)*0.5d0
	if(ipol.eq.2.and.nmo(2).ne.0) then
         xb=(evals(i+nbf)-efermi)/ssmear
         nb=util_erfc(xb)*0.5d0
        endif
         fzeroa=fzeroa+na
         if(ipol.eq.2) fzerob=fzerob+nb
      enddo
c      write(6,*) ' ef a b ',efermi,fzeroa,fzerob
         fzero=ntotel-(fzeroa+fzerob)*rhfuhf
c      write(6,*) ' ntotel  ',ntotel
c      write(6,*) ' fzero efer ',fzero,efermi
c      write(6,*) ' ==============='
      return
      end
c
      subroutine dft_focdm(focc,noc,geom,AO_bas_han,nbf,
     *     g_vecs,g_dens,toll)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
      integer noc
      integer AO_bas_han,nbf,geom,g_vecs,g_dens
      double precision focc(*),toll
c
      integer me,nproc,ltmpm,itmpm,g_tmp,i,nocold
      double precision zero
      parameter(zero=0.d0)
      logical debug
c
c     preliminaries
      me = ga_nodeid()
      nproc = ga_nnodes()
      nocold=noc
      debug = .false.
c
      if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &   call errquit('dftdensm: failed to alloc tmpm',0, MA_ERR)
      if (.not. ga_duplicate(g_dens, g_tmp, 'ga_temp'))
     &   call errquit('dft_densm: error creating ga',0,GA_ERR)
      call ga_zero(g_tmp)
      do i = 1, nbf
        if(focc(i).gt.toll) then
          noc=i
          if(mod(i-1+me,nproc).eq.0) then
            call get_col(g_vecs, nbf, i, DBL_MB(itmpm))
            call dscal(nbf, focc(i), DBL_MB(itmpm), 1)
            call put_col(g_tmp, nbf, i, DBL_MB(itmpm))
          endif
        endif
      enddo
      call ga_sync()
      call ga_dgemm('N', 'T', nbf, nbf, noc,1.d0,g_tmp,
     $   g_vecs, zero, g_dens)
      if (.not.ma_pop_stack(ltmpm))
     &   call errquit('dft_densm: cannot pop stack',0, MA_ERR)
      if (.not. ga_destroy(g_tmp)) call errquit
     &   ('dft_densm: could not destroy g_tmp', 0, GA_ERR)
c     
c     reset noc for future calls
c     
      noc=nocold
      return
      end
c
c     determine fractional occupations for fon calculations
c
      logical function dft_fon(rtdb,nmo_fon,ncore_fon,nel_fon,
     .     nbf,ntotel,focc,noc,ipol,me)
c
      implicit none
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer ipol
      integer nbf
      double precision focc(nbf,ipol)
      integer ntotel
      integer noc(2)
      integer me
c
      integer i
      double precision avg_fon, avg_fon2
      double precision nel_fon(4)
      integer nmo_fon(4), ncore_fon(2)
      integer ispin
      double precision ncheck
      logical debug
      integer nfilled
      logical do_avg_fon
      logical do_core_fon
c
      debug = .false.
      if (.not. rtdb_get(rtdb, 'dft:debugfon', mt_log, 1, debug)) 
     &     debug = .false.
c
c     do average fractional occupation by default
      do_avg_fon = .true.
      if (.not.rtdb_get(rtdb,'dft:avg_fon',mt_log,1,do_avg_fon))
     &     do_avg_fon = .true.
c
c     do average occupation starting with core orbitals
      do_core_fon = .false.
      if (.not.rtdb_get(rtdb,'dft:core_fon',mt_log,1,do_core_fon))
     &     do_core_fon = .false.
c
      if (me.eq.0 .and. debug) then
        write (luout,*) 'FON: ipol, noc, ntotel',ipol,noc(:),ntotel
        write (luout,*) 'do_avg_fon: ', do_avg_fon
      end if                    ! debug
c
      ncheck = 0d0
      do ispin = 1,ipol
        if (nmo_fon(ispin).lt.1) 
     &     call errquit('dft_densm:fon nmo_fon <1',1, INPUT_ERR)
        if (nel_fon(ispin).lt.0d0) 
     &     call errquit('dft_scf_so:fon nel_fon <0',1, INPUT_ERR)
c
        avg_fon = nel_fon(ispin)/dble(nmo_fon(ispin))  ! average occupation
        if (nmo_fon(ispin+2).lt.1.or.nel_fon(ispin+2).le.0d0) then
            avg_fon2 = 0.0d0
        else
            avg_fon2 = nel_fon(ispin+2)/dble(nmo_fon(ispin+2))  ! average occupation
        endif
        nfilled = ncore_fon(ispin)  ! number of filled orbitals for each spin
c
c       initialize
        do i = 1,nbf
          focc(i,ispin) = 0.0
        end do   
c
c       fill the molecular orbitals either starting from the core or valence states
        if (do_core_fon) then
c
c         partially filled molecular orbitals
          do i = 1,nmo_fon(ispin)   
            if (i > nbf) call errquit(
     &       'dft_densm:fon focc index exceeds nbf',1,INPUT_ERR)
            if (do_avg_fon) then
             focc(i,ispin) = avg_fon  ! assign average occupation
            else
             focc(i,ispin) = nel_fon(ispin) ! assign given fractional electron
            end if ! do_avg_fon
            ncheck = ncheck + focc(i,ispin)
            if (me.eq.0 .and. debug) 
     &        write(luout,*) i,ispin,focc(i,ispin)
          end do  ! nmo_fon
c
c         fully filled molecular orbitals        
          do i = 1,nfilled
            if (i+nmo_fon(ispin)>nbf) call errquit(
     &       'dft_densm:fon focc index exceeds nbf',1,INPUT_ERR)
            focc(i+nmo_fon(ispin),ispin) = 2d0/ipol
            ncheck = ncheck + focc(i+nmo_fon(ispin),ispin)
            if (me.eq.0 .and. debug) 
     &        write(luout,*) i,ispin,focc(i+nmo_fon(ispin),ispin)
          end do  ! nfilled

        else
c
c         fully filled molecular orbitals
          do i = 1,nfilled
            if (i>nbf) call errquit(
     &       'dft_densm:fon focc index exceeds nbf',1,INPUT_ERR)
            focc(i,ispin) = 2d0/ipol
            ncheck = ncheck + focc(i,ispin)
            if (me.eq.0 .and. debug) 
     &        write(luout,*) i,ispin,focc(i,ispin)
          end do  ! nfilled
c
c         partially filled molecular orbitals
          do i = 1,nmo_fon(ispin)   
            if (i+nfilled > nbf) call errquit(
     &       'dft_densm:fon focc index exceeds nbf',1,INPUT_ERR)
            if (do_avg_fon) then
             focc(i+nfilled,ispin) = avg_fon  ! assign average occupation
            else
             focc(i+nfilled,ispin) = nel_fon(ispin) ! assign given fractional electron
            end if ! do_avg_fon
            ncheck = ncheck + focc(i+nfilled,ispin)
            if (me.eq.0 .and. debug) 
     &        write(luout,*) i,ispin,focc(i+nfilled,ispin)
          end do  ! nmo_fon
c
c         second partially filled molecular orbitals
          do i = 1,nmo_fon(ispin+2)
            if (i+nfilled > nbf) call errquit(
     &       'dft_densm:fon focc index exceeds nbf',1,INPUT_ERR)
            if (do_avg_fon) then
             focc(i+nfilled+nmo_fon(ispin),ispin) = avg_fon2  ! assign average occupation
            else
             focc(i+nfilled+nmo_fon(ispin),ispin) = nel_fon(ispin+2) ! assign given fractional electron
            end if ! do_avg_fon
            ncheck = ncheck + focc(i+nfilled+nmo_fon(ispin),ispin)
            if (me.eq.0 .and. debug)
     &       write(luout,*) i,ispin,focc(i+nfilled+nmo_fon(ispin),ispin)
          end do  ! nmo_fon second partial
       end if  ! do_core_fon
      end do                    ! ispin
c
      if (me.eq.0 .and. debug) then
        write(luout,*) 'FON: focc:',focc(:,1)
        write(luout,*) 'ncheck:' , ncheck
        write(luout,*) 'ntotel:' , ntotel
      end if   
c
      if(abs(ncheck-dble(ntotel)).gt.1d-3 .and. me.eq.0) then
         dft_fon=.false.
      else
         dft_fon=.true.
      endif
c
      return
      end
c $Id$
