      subroutine grad2(d_ij,d_kl,d_ik,d_jl,d_il,d_jk,frc,g_dens,ldim,
     $                 basis,nproc,nat, scr, lscr, buf, lbuf, tol2e,
     $                 s34, oactive, oskel)
c$Id: grad2.F,v 1.12 1995-10-16 18:21:35 d3g681 Exp $

      implicit real*8 (a-h,o-z)

#include "sym.fh"
#include "mafdecls.h"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
#if defined(REPORT_LOOP)        
      logical reportit
      character*26 date_string
#endif
      logical oskel
      integer g_dens, basis
      double precision q4, q2
      integer geom

      logical status, oactive, oij, okl, oikjl, oschwarz12, osym

      dimension scr(lscr), buf(lbuf), frc ( 3, nat ), idatom(4),
     $          s34(nat)
      dimension d_ij(ldim,ldim),d_kl(ldim,ldim),d_ik(ldim,ldim),
     $          d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim)
      dimension oactive(nat)

      if (.not. bas_geom(basis, geom)) 
     $     call errquit('grad2: bad basis',0)

      iproc = ga_nodeid()

      q4 = 1.0d0                ! If not using symmetry

      ijatom = -1
      next = nxtval(nproc)
      do 95, iat3 = nat, 1, -1
#if defined(REPORT_LOOP)        
        if(ga_nodeid().eq.0) then
          reportit = mod(iat3,10).eq.0
          reportit = reportit.or.(iat3.eq.1)
          reportit = reportit.or.(iat3.eq.nat)
          if (reportit) then
            call util_date(date_string)
            write(6,*) 'grad2: atom =',iat3,' ',date_string
            call util_flush(6)
          endif
        endif
#endif
        do iat4 = 1, iat3
          s34(iat4) = schwarz_atom(iat3,iat4) ! Precompute
        enddo

        do 90, iat1 = nat, iat3, -1 ! Want canonical atom list for symmetry
          do 85, iat2 = iat1, 1, -1

            s12 = schwarz_atom(iat1,iat2)
            oschwarz12 =(s12*schwarz_max().ge.tol2e)
C
            if (oschwarz12 .and. oskel) then
               oschwarz12 = sym_atom_pair(geom, iat1, iat2, q2)
            endif
c
            if (oschwarz12) then
              ihi4 = iat3
              if ( iat1.eq.iat3 ) ihi4 = iat2
              do 80, iat4 = 1, ihi4

C                 check if all atoms are inactive
                  if ( .not.oactive(iat1) .and. .not.oactive(iat2) .and.
     $                 .not.oactive(iat3) .and. .not.oactive(iat4) )
     $                 goto 1020

                if ( s12 * s34(iat4) .lt. tol2e ) goto 1020

                if (oskel) then
                   osym = sym_atom_quartet(geom,
     $                  iat1, iat2, iat3, iat4, q4)
*                   write(6,11) iat1, iat2, iat3, iat4, q4, osym
* 11                format(1x,4i5,f6.1,2x,l1)
                   if (.not. osym) goto 1020
                endif

                status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2bfr',0)
                status = bas_ce2bfr(basis,iat2,iab2f,iab2l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2bfr',0)
                status = bas_ce2bfr(basis,iat3,iab3f,iab3l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2bfr',0)
                status = bas_ce2bfr(basis,iat4,iab4f,iab4l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2bfr',0)

                if (iab1f.le.0 .or. iab2f.le.0 .or. iab3f.le.0 .or.
     $               iab4f.le.0) then
c     
c     At least one center has no functions on it ... next atom
c
                   goto 1020
                endif

                status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2cnr',0)
                status = bas_ce2cnr(basis,iat2,iac2f,iac2l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2cnr',0)
                status = bas_ce2cnr(basis,iat3,iac3f,iac3l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2cnr',0)
                status = bas_ce2cnr(basis,iat4,iac4f,iac4l)
                if (.not.status) 
     $              call errquit('grad2: error with bas_ce2cnr',0)

              ijatom = ijatom + 1
              if ( next .eq. ijatom ) then

C            get density matrix blocks
             call get_dens ( d_ij, d_kl, d_ik, d_jl, d_il, d_jk, ldim,
     $           iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l,
     $           g_dens )

                oij = iat1.eq.iat2
                okl = iat3.eq.iat4
                oikjl = iat1.eq.iat3 .and. iat2.eq.iat4
                do 75, ish1 = iac1f, iac1l
                  if ( oij ) iac2l = ish1
                  do 70, ish2 = iac2f, iac2l
                    if (schwarz_shell(ish1,ish2)*schwarz_max().ge.tol2e) 
     $                  then
                      itop3 = iac3l
                      if ( iat1.eq.iat3 .and. iat2.eq.iat4 ) itop3 =ish1
                      do 65, ish3 = iac3f, itop3
                        itop4 = iac4l
                        if ( okl ) itop4 = ish3
                        if ( oikjl .and. ish1.eq.ish3 ) itop4 = ish2
                        do 60, ish4 = iac4f, itop4
                          if (schwarz_shell(ish1,ish2)*
     $                        schwarz_shell(ish3,ish4) .ge. tol2e) then
C---------------------------------------------------------------------

C     find contraction in atomic block

      status = bas_cn2bfr(basis,ish1,if1,il1)
      if1 = if1 - iab1f + 1
      il1 = il1 - iab1f + 1
      status = bas_cn2bfr(basis,ish2,if2,il2)
      if2 = if2 - iab2f + 1
      il2 = il2 - iab2f + 1
      status = bas_cn2bfr(basis,ish3,if3,il3)
      if3 = if3 - iab3f + 1
      il3 = il3 - iab3f + 1
      status = bas_cn2bfr(basis,ish4,if4,il4)
      if4 = if4 - iab4f + 1
      il4 = il4 - iab4f + 1

c Determine the number of ERIs in this block of integrals.

      nshbfp = ( il1 - if1 + 1 ) * ( il2 - if2 + 1 )
      nshbfq = ( il3 - if3 + 1 ) * ( il4 - if4 + 1 )
      Nint = nshbfp*nshbfq

C     2el. integral derivatives
      call intd_2e4c(basis, ish1, ish2, basis, ish3, ish4,
     &       lscr, scr, lbuf, buf, idatom)

      icount = 1
      do 38, i = 1, 4
        if ( idatom(i).ge.1 ) then
          do 35, icart = 1, 3
            dE = 0.D0
            do 55, ip1 = if1, il1
              do 50, ip2 = if2, il2
                do 45, ip3 = if3, il3
                  do 40, ip4 = if4, il4
                    den2 = 0.5D0  * ( d_ij(ip1,ip2) * d_kl(ip3,ip4) )
     $                   - 0.125D0 * (d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                             +  d_il(ip1,ip4) * d_jk(ip2,ip3) )
                    dE = dE + den2 * buf(icount)
                    icount = icount + 1

 40               continue
 45             continue
 50           continue
 55         continue
            if ( .not.oij .or. ish1.ne.ish2 ) dE = dE + dE
            if ( .not.okl .or. ish3.ne.ish4 ) dE = dE + dE
            if ( .not.oikjl .or. ish1.ne.ish3 .or. ish2.ne.ish4 )
     $           dE = dE + dE
*
*            if (icart.eq.3 .and. idatom(i).eq.1) then
*               write(6,17) iat1, iat2,iat3,iat4,dE, q4
* 17            format(1x,4i5,f10.6,f6.1)
*            endif
            frc(icart,idatom(i)) = frc(icart,idatom(i)) + dE * q4
 35       continue
        else
          icount = icount + 3 * Nint
        endif
 38   continue

C---------------------------------------------------------------------
                          endif   !schwarz
 60                     continue
 65                   continue
                    endif         !schwarz
 70               continue
 75             continue
                next = nxtval(nproc)
              endif ! next

 1020             continue
 80             continue
 83           continue
            endif
 85       continue
 90     continue
 95   continue
      next = nxtval(-nproc)
      return
      end


