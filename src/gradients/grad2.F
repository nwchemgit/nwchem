      subroutine grad2 (d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $                  d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2, 
     $                  frc, g_dens, ldim,
     $                  basis, nproc, nat, scr, lscr, buf, lbuf, tol2e,
     $                  s34, oactive, oskel, scftype)
c$Id: grad2.F,v 1.15 1995-11-13 16:56:39 d3g681 Exp $

      implicit none

#include "sym.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"

C-----------------------parameters---------------------------------
      integer ldim,          ! max dimension of density block
     $     g_dens,           ! array of global density handles
     $     basis, geom,      ! handles
     $     nproc, nat, lscr, lbuf

      dimension g_dens(3)

      double precision frc,  ! forces on atoms (local)
     $     scr, buf, s34, ! buffers and intermediates
     $     tol2e,
     $     d_ij, d_kl, d_ik, d_jl, d_il, d_jk,       ! local density
     $     d_ij2, d_kl2, d_ik2, d_jl2, d_il2, d_jk2  ! matrix blocks


      dimension scr(lscr), buf(lbuf), frc ( 3, nat ), s34(nat),
     $          d_ij(ldim,ldim),d_kl(ldim,ldim),d_ik(ldim,ldim),
     $          d_jl(ldim,ldim),d_il(ldim,ldim),d_jk(ldim,ldim),
     $          d_ij2(ldim,ldim), d_kl2(ldim,ldim),
     $          d_ik2(ldim,ldim), d_jl2(ldim,ldim),
     $          d_il2(ldim,ldim), d_jk2(ldim,ldim)


      logical oactive, oskel

      dimension oactive(nat)

      character*4 scftype

C-------------------------local variables-----------------------------
      integer i, ijatom, next, icart, idatom, icount,
     $     iat1, iat2, iat3, iat4, ihi4, ish1, ish2, ish3, ish4,
     $     iab1f, iab1l, iab2f, iab2l, iab3f, iab3l, iab4f, iab4l, 
     $     iac1f, iac1l, iac2f, iac2l, iac3f, iac3l, iac4f, iac4l,
     $     itop3, itop4, if1, il1, if2, il2, if3, il3, if4, il4,
     $     ip1, ip2, ip3, ip4, nshbfp, nshbfq, Nint

      integer nxtask, task_size
      external nxtask

      dimension idatom(4)

      double precision dE, den2, q4, q2, s12, half, quart, eigh, sixh

      logical status, oij, okl, oikjl, oschwarz12, osym, ouhf, orohf

      parameter ( half = 0.5D0,
     $            quart= 0.25D0,
     $            eigh = 0.125D0,
     $            sixh = 0.0625D0 )

      
#if defined(REPORT_LOOP)        
      logical reportit
      character*26 date_string
#endif

      if (.not. bas_geom(basis, geom)) 
     $     call errquit('grad2: bad basis',0)

      ouhf  = (scftype .eq. 'UHF')
      orohf = (scftype .eq. 'ROHF')

      q4 = 1.0d0                ! If not using symmetry

      call fock_2e_est_task(geom, oskel, tol2e, task_size)

      ijatom = -1
      next = nxtask(nproc,task_size)
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
                if ( ouhf ) then
                  call grad_get_dens_uhf ( d_ij, d_kl, d_ik, 
     $                               d_jl, d_il, d_jk,
     $                               d_ik2, d_jl2, d_il2, d_jk2, ldim, 
     $                               iab1f, iab1l, iab2f, iab2l, 
     $                               iab3f, iab3l, iab4f, iab4l,
     $                               g_dens )
                else if ( orohf ) then
                  call grad_get_dens_rohf ( d_ij, d_kl, d_ik, 
     $                               d_jl, d_il, d_jk,
     $                               d_ij2, d_kl2, d_ik2, 
     $                               d_jl2, d_il2, d_jk2, ldim, 
     $                               iab1f, iab1l, iab2f, iab2l, 
     $                               iab3f, iab3l, iab4f, iab4l,
     $                               g_dens )
                else
                  call grad_get_dens ( d_ij, d_kl, d_ik, 
     $                            d_jl, d_il, d_jk, ldim,
     $                            iab1f, iab1l, iab2f, iab2l, 
     $                            iab3f, iab3l, iab4f, iab4l,
     $                            g_dens(1) )
                endif
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
          if ( ouhf ) then
            do 135, icart = 1, 3
              dE = 0.D0
              do 155, ip1 = if1, il1
                do 150, ip2 = if2, il2
                  do 145, ip3 = if3, il3
                    do 140, ip4 = if4, il4
                      den2 = half * ( d_ij(ip1,ip2) * d_kl(ip3,ip4) )
     $                     - eigh  * (d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                             +  d_il(ip1,ip4) * d_jk(ip2,ip3)
     $                             +  d_ik2(ip1,ip3) * d_jl2(ip2,ip4)
     $                             +  d_il2(ip1,ip4) * d_jk2(ip2,ip3))
                      dE = dE + den2 * buf(icount)
                      icount = icount + 1

 140                continue
 145              continue
 150            continue
 155          continue
              if ( .not.oij .or. ish1.ne.ish2 ) dE = dE + dE
              if ( .not.okl .or. ish3.ne.ish4 ) dE = dE + dE
              if ( .not.oikjl .or. ish1.ne.ish3 .or. ish2.ne.ish4 )
     $             dE = dE + dE
              frc(icart,idatom(i)) = frc(icart,idatom(i)) + dE * q4
 135        continue
          else if ( orohf ) then
            do 235, icart = 1, 3
              dE = 0.D0
              do 255, ip1 = if1, il1
                do 250, ip2 = if2, il2
                  do 245, ip3 = if3, il3
                    do 240, ip4 = if4, il4
                      den2 = half * ( d_ij(ip1,ip2) * d_kl(ip3,ip4) )
     $                     + quart* (d_ij(ip1,ip2) * d_kl2(ip3,ip4)
     $                            +  d_ij2(ip1,ip2) * d_kl(ip3,ip4) )
     $                     + eigh * (d_ij2(ip1,ip2) * d_kl2(ip3,ip4) )
     $                     - eigh * (d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                            +  d_il(ip1,ip4) * d_jk(ip2,ip3) )
     $                     - sixh * (d_ik(ip1,ip3) * d_jl2(ip2,ip4)
     $                            +  d_il2(ip1,ip4) * d_jk(ip2,ip3)
     $                            +  d_ik2(ip1,ip3) * d_jl(ip2,ip4)
     $                            +  d_il(ip1,ip4) * d_jk2(ip2,ip3)
     $                            +  d_ik2(ip1,ip3) * d_jl2(ip2,ip4)
     $                            +  d_il2(ip1,ip4) * d_jk2(ip2,ip3) )
                      dE = dE + den2 * buf(icount)
                      icount = icount + 1

 240                continue
 245              continue
 250            continue
 255          continue
              if ( .not.oij .or. ish1.ne.ish2 ) dE = dE + dE
              if ( .not.okl .or. ish3.ne.ish4 ) dE = dE + dE
              if ( .not.oikjl .or. ish1.ne.ish3 .or. ish2.ne.ish4 )
     $             dE = dE + dE
              frc(icart,idatom(i)) = frc(icart,idatom(i)) + dE * q4
 235        continue
          else   ! RHF
            do 35, icart = 1, 3
              dE = 0.D0
              do 55, ip1 = if1, il1
                do 50, ip2 = if2, il2
                  do 45, ip3 = if3, il3
                    do 40, ip4 = if4, il4
                      den2 = half * ( d_ij(ip1,ip2) * d_kl(ip3,ip4) )
     $                     - eigh * (d_ik(ip1,ip3) * d_jl(ip2,ip4)
     $                            +  d_il(ip1,ip4) * d_jk(ip2,ip3) )
                      dE = dE + den2 * buf(icount)
                      icount = icount + 1

 40                 continue
 45               continue
 50             continue
 55           continue
              if ( .not.oij .or. ish1.ne.ish2 ) dE = dE + dE
              if ( .not.okl .or. ish3.ne.ish4 ) dE = dE + dE
              if ( .not.oikjl .or. ish1.ne.ish3 .or. ish2.ne.ish4 )
     $             dE = dE + dE
              frc(icart,idatom(i)) = frc(icart,idatom(i)) + dE * q4
 35         continue
          endif   ! rhf
        else  ! idatom
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
                next = nxtask(nproc,task_size)
              endif ! next

 1020             continue
 80             continue
 83           continue
            endif
 85       continue
 90     continue
 95   continue
      next = nxtask(-nproc,task_size)
      return
      end


