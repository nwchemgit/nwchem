      subroutine dftg_grad2(d_ija,d_kla,frc,g_dens,
     &                 ldim,
     $                 basis,nproc,nat, scr, lscr, buf, lbuf,
     &                 tol2e, oactive, nactive )
c****
c**** exchange density matrices removed from call list for DFT
c**** coulomb-only term.
c****

      implicit real*8 (a-h,o-z)

#include "mafdecls.h"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "tcgmsg.fh"

      integer g_dens(2), basis, nactive

      logical status, oij, okl,oikjl, oactive(nat)

      dimension scr(lscr), buf(lbuf), frc ( 3, nat ), idatom(4)
      dimension d_ija(ldim,ldim),d_kla(ldim,ldim)

      iproc = ga_nodeid()

      ijatom = -1
      next = nxtval(nproc)
      do 95, iat3 = nat, 1, -1
        do 90, iat1 = iat3, 1, -1
          do 85, iat2 = iat1, 1, -1
            if (schwarz_atom(iat1,iat2)*schwarz_max() .ge. tol2e) then
              ihi4 = iat3
              if ( iat1.eq.iat3 ) ihi4 = iat2
              do 80, iat4 = ihi4, 1, -1

                if ( schwarz_atom(iat1,iat2) * schwarz_atom(iat3,iat4) 
     $               .ge. tol2e ) then

              ijatom = ijatom + 1
C             if (mod(ijatom, nproc) .eq. iproc) then
              if ( next .eq. ijatom ) then

c****           include only terms with at least one active center
                if ( .not.oactive(iat1) .and. .not.oactive(iat2)
     &               .and.
     &               .not.oactive(iat3) .and. .not.oactive(iat4) ) 
     &               then 
c                 All centers are inactive; next atom
                  goto 1010
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

                if (iab1f.le.0 .or. iab2f.le.0 .or. iab3f.le.0 .or.
     $               iab4f.le.0) then
c     
c     At least one center has no functions on it ... next atom
c
                   goto 1010
                endif

C               get density matrix blocks
c****           
c****           modified for coulomb-only and UHF/RHF.  For exact exchange
c****           need to use get_dens calls (for UHF, 2 calls)
c****           
                call dftg_get_cden ( d_ija, d_kla,  ldim,
     $               iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, 
     $               iab4f, iab4l, g_dens(1) )
c****

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
      do i = 1, 4
        if ( idatom(i).ge.1 ) then
          do icart = 1, 3
            dE = 0.D0
            do ip1 = if1, il1
              do ip2 = if2, il2
                do ip3 = if3, il3
                  do ip4 = if4, il4
c****               
c****               hacked for DFT; exchange contributions deleted
c****               
c****               den2 = 0.5D0  * ( d_ij(ip1,ip2) * d_kl(ip3,ip4) )
c****               $             - 0.125D0 * (d_ik(ip1,ip3) * d_jl(ip2,ip4)
c****               $                       +  d_il(ip1,ip4) * d_jk(ip2,ip3) )

                    den2 = 0.5D0 *
     &                   (d_ija(ip1,ip2) * d_kla(ip3,ip4))

                    dE = dE + den2 * buf(icount)
                    icount = icount + 1

                  enddo 
                enddo
              enddo
            enddo

            if ( .not.oij .or. ish1.ne.ish2 ) dE = dE + dE
            if ( .not.okl .or. ish3.ne.ish4 ) dE = dE + dE
            if ( .not.oikjl .or. ish1.ne.ish3 .or. ish2.ne.ish4 )
     $           dE = dE + dE
            frc(icart,idatom(i)) = frc(icart,idatom(i)) + dE
          enddo
        else
          icount = icount + 3 * Nint
        endif
      enddo 
C---------------------------------------------------------------------
                          endif   !schwarz
 60                     continue
 65                   continue
                    endif         !schwarz
 70               continue
 75             continue
 1010           next = nxtval(nproc)
              endif

                  endif
 80             continue
 83           continue
            endif
 85       continue
 90     continue
 95   continue
      next = nxtval(-nproc)
      return
      end
