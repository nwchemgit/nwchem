#ifndef SECOND_DERIV
      Subroutine dftg_cdfit (rtdb, geom, AO_bas_han, CD_bas_han,
     &                       nbf_ao, nbf_cd, nat, tol2e, scr, 
     &                       lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                       max_sh_bfcd,
     &                       iga_dens, force, CD_coef, oskel)
#else
      Subroutine dfth_cdfit (rtdb, geom, AO_bas_han, CD_bas_han,
     &                       nbf_ao, nbf_cd, nat, tol2e, scr, 
     &                       lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                       max_sh_bfcd,
     &                       iga_dens, hess, CD_coef, oskel)
#endif
*
* $Id: dftg_cdfit.F,v 1.14 1999-02-15 18:27:42 bjohnson Exp $
*
      implicit none
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "util.fh"
#include "sym.fh"
c
      integer rtdb, geom, AO_bas_han, CD_bas_han
      integer nbf_ao, nbf_cd, nat, max_sh_bf,max_sh_bfcd
      integer iga_dens, lbuf, lscr
c
c     scratch arrays for calls to integrals package
c
      double precision scr(lscr), buf(lbuf)
      integer idatom(4)
c
c     scratch arrays for dgem
c
      double precision vec(max_sh_bfcd*max_sh_bfcd),
     &                 PPP(max_sh_bf,max_sh_bf)
c
c-->  Charge Density Expansion Coefficients
c     
      double precision CD_coef(nbf_cd)
#ifndef SECOND_DERIV
c     
c     force vector
c     
      double precision force(3,nat)
#else
c     
c     hessian matrix
c     
      double precision hess(3,nat,3,nat)
#endif
c
      double precision tol2e
      logical oskel
c
c     local declarations
c
      logical doesit
      integer nproc,me,icount,igran
      integer nshells_ao,nshells_cd,ishp,ishq
      integer ifirstq, ilastq, ifirstp, ilastp,nshbfq,nshbfp,Nao2
      integer mu,nu,Nintegrals,ishd,ifirstd,ilastd,nshbfd
      integer iat,icount2,ifp,i,next,nxtask,icart,ist
      double precision PPP_max,dE,fac,schwarz_shell,q1
c
      character*(nw_max_path_len) cdfit
c
      external nxtask,schwarz_shell
      nproc  = ga_nnodes()
      me = ga_nodeid()
c
c     read CD_coeff vector written by fitcd
c
      if (me.eq.0)then
         call util_file_name('cdfit', .true., .false., cdfit)
         inquire(file=cdfit,exist=doesit)
         if (.not. doesit)call errquit
#ifndef SECOND_DERIV
     &        ('dftg_cdfit: file cdfit does not exist',0)
#else
     &        ('dfth_cdfit: file cdfit does not exist',0)
#endif
         open(unit=79,file=cdfit,access='sequential',
     &        form='unformatted',status='old')
         read(79)(CD_coef(i),i=1,nbf_cd)
         close(79)
      endif
c
      call ga_brdcst(Msg_dftg_cdf,CD_coef,mdtob(nbf_cd),0)
c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
      if (.not. bas_numcont(AO_bas_han,nshells_ao))then
#ifndef SECOND_DERIV
         call errquit('Exiting in dftg_cdfit.',2)
#else
         call errquit('Exiting in dfth_cdfit.',2)
#endif
      endif
c      
      if (.not. bas_numcont(CD_bas_han,nshells_cd))then
#ifndef SECOND_DERIV
         call errquit('Exiting in dftg_cdfit.',3)
#else
         call errquit('Exiting in dfth_cdfit.',3)
#endif
      endif
c
c     2el. 3-c integral derivatives
c
c     increased granularity
c
      igran = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      igran = max(1,igran)
c
      next = nxtask(nproc,igran) + 1
c
      do 205 ishp = 1, nshells_ao
         if (ishp.eq.next)then
c
            next = nxtask(nproc,igran) + 1
c
            if (.not. bas_cn2bfr(AO_bas_han,ishp,ifirstp,ilastp))
#ifndef SECOND_DERIV
     &         call errquit('Exiting in dftg_cdfit.',6)
#else
     &         call errquit('Exiting in dfth_cdfit.',6)
#endif
            nshbfp = ilastp - ifirstp + 1
c          
            do ishq = 1, ishp
               if (.not. bas_cn2bfr(AO_bas_han,ishq,ifirstq,ilastq))
#ifndef SECOND_DERIV
     &            call errquit('Exiting in dftg_cdfit.',7)
#else
     &            call errquit('Exiting in dfth_cdfit.',7)
#endif
               nshbfq = ilastq - ifirstq + 1
               Nao2 = nshbfp*nshbfq
c
               fac = 2d0
               if (ishp.eq.ishq) fac = 1d0
c
c              get shell block of DM
c
               call ga_get (iga_dens, ifirstq, ilastq, ifirstp, ilastp,
     &                      PPP, max_sh_bf)
               
c
c              Schwarz screening on product MAX(ABS(DM(i,j)))*schwarz(i,j)
c
               PPP_max = 0.d0
               do nu = 1, nshbfp
                  do mu = 1, nshbfq
                     PPP_max = max(PPP_max,abs(PPP(mu,nu)))
                  enddo
               enddo
               if (PPP_max*schwarz_shell(ishp,ishq).gt.tol2e)then
                  do ishd = 1, nshells_cd
                     if (.not. bas_cn2bfr(CD_bas_han, ishd, 
     &                  ifirstd, ilastd))
#ifndef SECOND_DERIV
     &                  call errquit('Exiting in dftg_cdfit.',8)
#else
     &                  call errquit('Exiting in dfth_cdfit.',8)
#endif
                     nshbfd = ilastd - ifirstd + 1
                     Nintegrals = nshbfp*nshbfq*nshbfd
c
c                    check for use of symmetry
c
                     if (oskel) then
                        q1 = 0.d0
                        if (sym_shell(cd_bas_han, ishd, q1))then
#ifndef SECOND_DERIV
                           call intd_2e3c(CD_bas_han, ishd, 
     &                                    AO_bas_han, ishp, ishq,
     &                                    lscr, scr, lbuf, buf, idatom)
#else
               write(*,*)'Call intdd_2e3c here!'
               call dfill(lbuf, 0.0d0, buf, 1)
               call ifill(4, 0, idatom, 1)
c                           call intdd_2e3c(CD_bas_han, ishd, 
c     &                                     AO_bas_han, ishp, ishq,
c     &                                     lscr, scr, lbuf, buf, idatom)
#endif
                        endif
                     else
                        q1 = 1.0d0
#ifndef SECOND_DERIV
                        call intd_2e3c(CD_bas_han, ishd, 
     &                                 AO_bas_han, ishp, ishq,
     &                                 lscr, scr, lbuf, buf, idatom)
#else
               write(*,*)'Call intdd_2e3c here!'
               call dfill(lbuf, 0.0d0, buf, 1)
               call ifill(4, 0, idatom, 1)
c                        call intdd_2e3c(CD_bas_han, ishd, 
c     &                                  AO_bas_han, ishp, ishq,
c     &                                  lscr, scr, lbuf, buf, idatom)
#endif
c
                     endif                   
c
                     if(q1.eq.0.d0)goto 204
#ifndef SECOND_DERIV
                     icount = 1
                     do iat = 1, 4
c
                        if (idatom(iat).gt.0)then
c
                           do icart = 1, 3
c
                              call dgemv('N', Nao2, nshbfd, 1.D0, 
     &                                   buf(icount), Nao2,
     &                                   CD_coef(ifirstd), 1, 
     &                                   0.D0, vec, 1)
c
                              icount = icount + Nintegrals
                              dE = 0.d0
                              icount2 = 1
c
                              do ifp = 1, nshbfp
                                 dE = dE + ddot(nshbfq, PPP(1,ifp), 
     &                                          1, vec(icount2), 1)
                                 icount2 = icount2 + nshbfq
                              enddo
c
                              force(icart,idatom(iat)) = 
     &                        force(icart,idatom(iat)) + dE * fac * q1
c
                           enddo
                        else
                           icount = icount + 3*Nintegrals
                        endif
                     enddo
#else
                     write(*,*)'Are dimensions right on ints and P?'
                     call d2int_2e3c_block_jhess(CD_coef(ifirstd), PPP,
     &                    max_sh_bf, hess, nat, buf, fac*q1, idatom,
     &                    nshbfd, nshbfp, nshbfq)
#endif
  204                continue
                  enddo
               endif
            enddo
         endif
  205 continue
c
      next = nxtask(-nproc,igran)
c
c     2el. 2-c integral derivatives
c
      igran = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      igran = max(1,igran)
c
      next = nxtask(nproc,igran) + 1
c
      do 305 ishp = 1, nshells_cd
c
         if (ishp.eq.next)then
c
            next = nxtask(nproc,igran) + 1
c
            if (.not. bas_cn2bfr(CD_bas_han,ishp,ifirstp,ilastp))
#ifndef SECOND_DERIV
     &         call errquit('Exiting in dftg_cdfit.',9)
#else
     &         call errquit('Exiting in dfth_cdfit.',9)
#endif
               nshbfp = ilastp - ifirstp + 1
c
            do ishq = 1, ishp
c
               if (.not. bas_cn2bfr(CD_bas_han,ishq,ifirstq,ilastq))
#ifndef SECOND_DERIV
     &            call errquit('Exiting in dftg_cdfit.',10)
#else
     &            call errquit('Exiting in dfth_cdfit.',10)
#endif
               nshbfq = ilastq - ifirstq + 1
               Nintegrals = nshbfp*nshbfq
c
#ifndef SECOND_DERIV
               call intd_2e2c(CD_bas_han, ishp, CD_bas_han, ishq, 
     &              lscr, scr, lbuf, buf, idatom)
#else
               write(*,*)'Call intdd_2e2c here!'
c               call intdd_2e2c(CD_bas_han, ishp, CD_bas_han, ishq, 
c     &              lscr, scr, lbuf, buf, idatom)
#endif
c
               fac = -1.d0
               if (ishp.eq.ishq)fac=-.5d0
c
#ifndef SECOND_DERIV
               icount = 1
               do iat = 1, 2
                  if (idatom(iat).gt.0)then
                     ist=1
                     do icart = 1, 3
                        call dgemv('n', nshbfq, nshbfp,
     &                             1.d0, buf(icount), nshbfq,
     &                             CD_coef(ifirstp), 1,
     &                             0.d0, vec(ist), 1)
                        ist = ist + nshbfq
                        icount = icount + Nintegrals
                     enddo
                     call dgemv('t', nshbfq, 3,
     &                          fac, vec, nshbfq,
     &                          CD_coef(ifirstq), 1,
     &                          1.d0, force(1,idatom(iat)), 1)
                  else
                     icount = icount + 3*Nintegrals
                  endif
               enddo
#else
               write(*,*)'Are dimensions right on p and q?'
               call d2int_2e2c_block_jhess(CD_coef(ifirstp),
     &              CD_coef(ifirstq), hess, nat, buf, fac, idatom,
     &              nshbfp, nshbfq)
#endif
            enddo
         endif
  305 continue
c
      next=nxtask(-nproc,igran)
c
      return 
      end 
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "dftg_cdfit.F"
#endif
