      subroutine dft_densm(g_dens,g_vecs,nbf,nmo,ipol,
     ,     geom,AO_bas_han,noc,ntotel,
     ,	evals,focc,ssmear,tdots,
     ,     iter,search0,
     .     fon, nel_fon,nmo_fon,
     .     spinset)
C$Id: dft_densm.F,v 1.31 2001-12-13 02:19:01 edo Exp $
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"

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
      integer nel_fon(2),nmo_fon(2)

c
      integer i,isp,g_tmp
      integer me, nproc
      double precision zero,toll,efermi(2)
      double precision cksum,ntot(2)
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ltmpm,itmpm,ipoint
      double precision sqrtpi,kbau,eps,rhfuhf
      double precision one,x,xb
      double precision util_erfc
      external util_erfc
      logical dft_checkdg
      external dft_checkdg
      parameter (zero=0.d0,toll=1.d-9,one=1.d0,
     ,   kbau=1.d0,eps=1.d-4)
c
      sqrtpi=sqrt(acos(-1d0))
      me=ga_nodeid()
      nproc=ga_nnodes()
      rhfuhf=2d0
      if(ipol.eq.2) rhfuhf=1d0
      if(ssmear.lt.toll.or.iter.lt.-1) then
c
c     check degener
c     
         if (fon) then
            if(dft_checkdg(nmo_fon,nel_fon,
     .          nbf,ntotel,focc,noc,ipol,me)) then
                do isp=1,ipol
                   call dft_focdm(focc(1+(isp-1)*nbf),noc(isp),geom,
     ,                 AO_bas_han,nbf,
     *                 g_vecs(isp),g_dens(isp),toll)
                enddo
                return
             endif
          endif  
cedo         if (check_deg.and.ipol.eq.2) then
         do isp=1,ipol
            call ga_dgemm('n', 't', nbf, nbf, 
     $           noc(isp), 2d0/dble(ipol), g_vecs(isp),
     $           g_vecs(isp), zero, g_dens(isp))
         enddo
      else	
         if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &        call errquit('dftdensm: failed to alloc tmpm',0)
c
c     initialize ef
c
           efermi(1)=0.d0
           efermi(2)=0.d0
           g_tmp = 
     =          ga_create_atom_blocked(geom, AO_bas_han, 'ga_temp')
           call ga_zero(g_tmp)
         if(search0)then
           call dfill(nbf*ipol, 0.d0, focc, 1)
           if(spinset) then
              nmo(2)=0
              call dft_zero(2,nbf,nmo,noc(1),efermi(1),evals,
     ,             ssmear,toll)
              call dft_zero(2,nbf,nmo,noc(2),efermi(2),evals(nbf+1),
     .             ssmear,toll)
              nmo(2)=nmo(1)
           else
              call dft_zero(ipol,nbf,nmo,ntotel,efermi(1),evals,ssmear,
     .             toll)
              efermi(2)=efermi(1)
           endif
         endif
         do isp=1,ipol
            ntot(isp)=0.d0
            ipoint=(isp-1)*nbf
            do i = 1, nmo(isp)
               x=(evals(ipoint+i)-efermi(isp))/ssmear
              if(search0)
     *             focc(i+ipoint)= util_erfc(x)*0.5d0*rhfuhf
               if(focc(i+ipoint).gt.toll) then
                  noc(isp)=i
               else
                  focc(i+ipoint)=0.d0
               endif
            enddo
            call ga_zero(g_tmp)
            do i=1+me,nmo(isp),nproc
               if(focc(i+ipoint).gt.toll) then
                  call get_col(g_vecs(isp), nbf, i, DBL_MB(itmpm))
                  call dscal(nbf, focc(i+ipoint), DBL_MB(itmpm), 1)
                  call put_col(g_tmp, nbf, i, DBL_MB(itmpm))
               endif
            enddo
            do i=1,nmo(isp)
               ntot(isp)=ntot(isp)+focc(i+ipoint)
            enddo
            call ga_sync()
!            call ga_dgemm('n', 't', nbf, nmo(isp), noc(isp),1.d0,g_tmp,
            call ga_dgemm('n', 't', nbf, nbf, noc(isp),1.d0,g_tmp,
     $           g_vecs(isp), zero, g_dens(isp))
         enddo                  ! end of loop over spins
         if (.not.ma_pop_stack(ltmpm))
     &        call errquit('dft_densm: cannot pop stack',0)
           if (.not. ga_destroy(g_tmp)) call errquit
     &          ('dft_densm: could not destroy g_tmp', 0)
c     
c     compute entropy correction to total e
c     
         if(.not.fon) then
         cksum=0.d0
         do i=1,nmo(1)
            x=(evals(i)-efermi(1))/ssmear
            cksum=cksum+exp(-x*x)/sqrtpi
            if(ipol.eq.2) then
               xb=(evals(i+nbf)-efermi(2))/ssmear
               cksum=cksum+exp(-xb*xb)/sqrtpi
            endif
         enddo
         tdots=tdots-cksum*ssmear*0.5d0*rhfuhf
         if(ga_nodeid().eq.0) then 
            write(LuOut,1688) ssmear,
     ,             efermi(1),tdots,(i,ntot(i),i=1,ipol)
 1688       format(' sigma= ',d7.2,
     ,             ' eFermi= ',f8.5,' T*S= ',f9.6,
     ,           2(' nel(',i1,')=',f8.2)/) 
C            write(LuOut,1723) (focc(i),i=1,nbf*ipol)
C 1723       format(' Foccs ',(
         endif
         endif
      endif
      if(spinset) then
         noc(1)=nint(ntot(1))
         if(ipol.eq.2) noc(2)=nint(ntot(2))
      endif
      return
      end
      subroutine dft_zero(ipol,nbf,nmo,
     ,     ntotel,efermi,evals,ssmear,toll)
      implicit none
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer ipol
      integer nbf,nmo(2)
      double precision evals(nbf*ipol) 
      double precision ssmear
      double precision toll
      integer ntotel
      double precision efermi
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
         efermi=evals(ntotel/2)
      else
c
c     open shell
c
         efermi=0.5d0*(evals(ntotel/2)+evals(nbf+ntotel/2))
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
       call errquit('dft_zero: extrema:too many attempts',mx_attem)     
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
      double precision function fzero(
     *     ipol,efermi,evals,nbf,nmo,ntotel,ssmear)
      implicit none
      integer ipol
      integer nbf,nmo(2)
      double precision efermi
      double precision evals(nbf*ipol)
      integer ntotel 
      double precision ssmear
#include "msgids.fh"
c
      integer i
      double precision util_erfc
      external util_erfc
      double precision one,x,xb,rhfuhf,
     , fzeroa ,fzerob ,na,nb
      parameter(one=1.d0)
      rhfuhf=2d0/ipol
      fzeroa=0d0
      fzerob=0d0
      nb=0d0
      do i=1,nmo(1)
         x=(evals(i)-efermi)/ssmear
         xb=(evals(i+nbf)-efermi)/ssmear
         na=util_erfc(x)*0.5d0
	if(ipol.eq.2.and.nmo(2).ne.0) 
     *    nb=util_erfc(xb)*0.5d0
         fzeroa=fzeroa+na
         if(ipol.eq.2) fzerob=fzerob+nb
      enddo
C      write(0,*) ' ef a b ',efermi,fzeroa,fzerob
C      write(0,*) 'abs  a b ',abs(noc(1)-fzeroa),abs(noc(2)-fzerob)
         fzero=ntotel-(fzeroa+fzerob)*rhfuhf
C      write(0,*) ' fzero efer ',fzero,efermi
C      write(0,*) ' ==============='
      return
      end
      subroutine dft_focdm(focc,noc,geom,AO_bas_han,nbf,
     *     g_vecs,g_dens,toll)
      implicit none
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
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      parameter(zero=0.d0)
         me = ga_nodeid()
         nproc = ga_nnodes()
         nocold=noc
         if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &        call errquit('dftdensm: failed to alloc tmpm',0)
         g_tmp = 
     =        ga_create_atom_blocked(geom, AO_bas_han, 'ga_temp')
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
            call ga_dgemm('n', 't', nbf, nbf, noc,1.d0,g_tmp,
     $           g_vecs, zero, g_dens)
         if (.not.ma_pop_stack(ltmpm))
     &        call errquit('dft_densm: cannot pop stack',0)
         if (.not. ga_destroy(g_tmp)) call errquit
     &        ('dft_densm: could not destroy g_tmp', 0)
c
c        reset noc for future calls
c
       noc=nocold
c      call dfill(nbf, 0.d0, focc, 1)
c      call dfill(noc, rhfuhf, focc, 1)
         return
         end
      logical function dft_checkdg(nmo_fon,nel_fon,
     .     nbf,ntotel,focc,noc,ipol,
     , me)
      implicit none
#include "global.fh"
#include "stdio.fh"
      integer ipol
      integer nbf
      double precision focc(nbf,ipol)
      integer ntotel
      integer noc(2)
      integer me
c
      integer i
      double precision avg_fon
      integer nmo_fon(2),nel_fon(2),ispin
      double precision ncheck
c
c  assignments      
c
      ncheck=0d0
      do ispin=1,ipol
         avg_fon = dble(nel_fon(ispin))/dble(nmo_fon(ispin)) 
         do i=1,noc(ispin)-nel_fon(ispin)
            focc(i,ispin) = 2d0/ipol
            ncheck=ncheck+focc(i,ispin)
         enddo
         do i = noc(ispin)-nel_fon(ispin)+1,
     ,        noc(ispin)-nel_fon(ispin)+nmo_fon(ispin)
            focc(i,ispin) = avg_fon*(2d0/ipol)
            ncheck=ncheck+focc(i,ispin)
         enddo
      enddo
      if(abs(ncheck-dble(ntotel)).gt.1d-3) then
         write(0,*) ' lost electrons ',ncheck,' vs ',ntotel
         dft_checkdg=.false.
      else
         dft_checkdg=.true.
      endif
      if(me.eq.0) write(luout,*)  ' FON applied '
      return
      end
