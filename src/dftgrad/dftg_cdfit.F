      Subroutine dftg_cdfit(rtdb, geom, AO_bas_han, CD_bas_han,
     &                      nbf_ao, nbf_cd, nat, tol2e, scr, 
     &                      lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                      iga_dens, force, CD_coef,
     &                      oskel)
c
c     Turned into wrapper to call generalized function
c
c     BGJ - 3/99
c    
*
* $Id: dftg_cdfit.F,v 1.18 2000-02-09 22:27:45 edo Exp $
*
      integer rtdb, geom, AO_bas_han, CD_bas_han
      integer nbf_ao, nbf_cd, nat, max_sh_bf
      integer iga_dens, lbuf, lscr
      double precision scr(lscr), buf(lbuf)
      double precision vec(*),
     &                 PPP(max_sh_bf,max_sh_bf)
      double precision CD_coef(nbf_cd)
      double precision force(3,nat)
      double precision tol2e
      logical oskel
c
      call dftg_cdfit_gen(rtdb, geom, AO_bas_han, CD_bas_han,
     &                    nbf_ao, nbf_cd, nat, tol2e, scr, 
     &                    lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                    iga_dens, force, 0d0, 0, 0d0,
     &                    CD_coef, oskel, 1)
c
      return
      end
c
c     Generalized from dftg_cdfit to calculate fitted J derivative
c     matrices and hessian as well as fitted J forces
c
c     BGJ - 3/99
c
      Subroutine dftg_cdfit_gen(rtdb, geom, AO_bas_han, CD_bas_han,
     &                          nbf_ao, nbf_cd, nat, tol2e, scr, 
     &                          lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                          iga_dens, force, hess,
     &                          g_dj, db, CD_coef, oskel_in, calc_type)
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
      integer nbf_ao, nbf_cd, nat, max_sh_bf
      integer iga_dens, lbuf, lscr, calc_type
c
c     scratch arrays for calls to integrals package
c
      double precision scr(lscr), buf(lbuf)
      integer idatom(4)
c
c     scratch arrays for dgem
c
      double precision vec(*),
     &                 PPP(max_sh_bf,max_sh_bf)
c
c-->  Charge Density Expansion Coefficients
c     
      double precision CD_coef(nbf_cd)
c     
c     force vector
c     
      double precision force(3,nat)
c     
c     hessian matrix
c     
      double precision hess(3,nat,3,nat)
c
c     Handles for derivative J matrices
c
      integer g_dj(3,nat)
c
c     Right side vectors of fit required to calculate J matrix derivatives
c
      double precision db(nbf_cd,3,nat)
c
      double precision tol2e
      logical oskel_in
c
c     local declarations
c
      logical doesit, oskel
      integer nproc,me,icount,igran
      integer nshells_ao,nshells_cd,ishp,ishq
      integer ifirstq, ilastq, ifirstp, ilastp,nshbfq,nshbfp,Nao2
      integer mu,nu,Nintegrals,ishd,ifirstd,ilastd,nshbfd
      integer iat,icount2,ifp,i,next,nxtask,icart,ist,l_dfkl,k_dfkl
      integer iatom,ishatf,ishatl
      double precision PPP_max,dE,fac,schwarz_shell,q1
      logical do_force, do_cpks_r, do_hess
c
      character*(nw_max_path_len) cdfit
c
      external nxtask,schwarz_shell
      nproc  = ga_nnodes()
      me = ga_nodeid()
c
      do_force = calc_type .eq. 1
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_force .or. do_cpks_r .or. do_hess)) then
         call errquit('dftg_cdfit_gen: illegal calculation type',1)
      endif
c
c     Can't use symmetry as it stands for CPKS RHS
c
      if (do_cpks_r) then
         oskel = .false.
      else
         oskel = oskel_in
      endif
c
c     Need some additional scratch and initialization for CPKS RHS
c
      if (do_cpks_r) then
         if (.not.MA_push_get(MT_DBL, max_sh_bf*max_sh_bf, 'dfkl',
     &        l_dfkl, k_dfkl)) call
     &        errquit('dftg_cdfit_gen: cannot allocate dfkl',1)
         call dfill(nbf_cd*3*nat, 0.0d0, db, 1)
      endif
c
c     read CD_coeff vector written by fitcd
c
      if (me.eq.0)then
         call util_file_name('cdfit', .true., .false., cdfit)
         inquire(file=cdfit,exist=doesit)
         if (.not. doesit)call errquit
     &        ('dftg_cdfit_gen: file cdfit does not exist',0)
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
         call errquit('Exiting in dftg_cdfit_gen.',2)
      endif
c      
      if (.not. bas_numcont(CD_bas_han,nshells_cd))then
         call errquit('Exiting in dftg_cdfit_gen.',3)
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
     &         call errquit('Exiting in dftg_cdfit_gen.',6)
            nshbfp = ilastp - ifirstp + 1
c          
            do ishq = 1, ishp
               if (.not. bas_cn2bfr(AO_bas_han,ishq,ifirstq,ilastq))
     &            call errquit('Exiting in dftg_cdfit_gen.',7)
               nshbfq = ilastq - ifirstq + 1
               Nao2 = nshbfp*nshbfq
c
               if (do_force.or.do_hess) then
                  fac = 2d0
                  if (ishp.eq.ishq) fac = 1d0
               else
                  fac = 1d0
               endif
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
     &                  call errquit('Exiting in dftg_cdfit_gen.',8)
                     nshbfd = ilastd - ifirstd + 1
                     Nintegrals = nshbfp*nshbfq*nshbfd
c
c                    check for use of symmetry
c
                     if (oskel) then
                        q1 = 0.d0
                        if (sym_shell(cd_bas_han, ishd, q1))then
                           if (do_force.or.do_cpks_r) then
                              call intd_2e3c(CD_bas_han, ishd, 
     &                             AO_bas_han, ishp, ishq,
     &                             lscr, scr, lbuf, buf, idatom)
                           else if (do_hess) then
                              call intdd_2e3c(CD_bas_han, ishd, 
     &                             AO_bas_han, ishp, ishq,
     &                             lscr, scr, lbuf, buf, idatom)
                           endif
                        endif
                     else
                        q1 = 1.0d0
                        if (do_force.or.do_cpks_r) then
                           call intd_2e3c(CD_bas_han, ishd, 
     &                          AO_bas_han, ishp, ishq,
     &                          lscr, scr, lbuf, buf, idatom)
                        else if (do_hess) then
c                           write(*,*)'Calling intdd_2e3c! no symm',
c     '                         ishd,ishp,ishq
                           call intdd_2e3c(CD_bas_han, ishd, 
     &                          AO_bas_han, ishp, ishq,
     &                          lscr, scr, lbuf, buf, idatom)
                        endif
c
                     endif                   
c
                     if(q1.eq.0.d0)goto 204
c
                     if (do_force) then
c
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
c
                     else if (do_hess) then
c
                        write(*,*)'Are dimensions right on ints and P?'
                        call d2int_2e3c_block_jhess(CD_coef(ifirstd),
     &                       PPP, max_sh_bf, hess, nat, buf, fac*q1,
     &                       idatom, nshbfd, nshbfp, nshbfq)
c
                     else if (do_cpks_r) then
c
                        write(*,*)'Are dimensions right on ints and P?'
                        call dint_2e3c_block_jdmat(CD_coef(ifirstd),
     &                       PPP, nbf_cd, max_sh_bf, g_dj,
     &                       db(ifirstd,1,1), dbl_mb(k_dfkl), buf,
     &                       ishp .eq. ishq, idatom, nshbfd,
     &                       ifirstp, ilastp, ifirstq, ilastq)
c
                     endif
  204                continue
                  enddo
               endif
            enddo
         endif
  205 continue
c
      next = nxtask(-nproc,igran)
c
      if (do_cpks_r) then
         if (.not.MA_pop_stack(l_dfkl))
     &        call errquit('dftg_cdfit_gen: cannot pop dfkl',1)
      endif
      if (do_cpks_r) then
         write(*,*)'DB after 3c part:'
         call output(db, 1, nbf_cd, 1, 3*nat, nbf_cd, 3*nat, 1)
c!!!         call dfill(nbf_cd*3*nat, 0.0d0, db, 1)
      endif
c
c     2el. 2-c integral derivatives
c
      igran = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      igran = max(1,igran)
c
      next = nxtask(nproc,igran) + 1
c
c     Since 2e2c integral derivatives are zero unless the centers
c     involved are different, the loops are set up to run only over
c     pairs of shells on different centers, saving a small amount
c     of work
c
      do 304 iatom = 1, nat
         if (.not. bas_ce2cnr(CD_bas_han,iatom,ishatf,ishatl))
     &        call errquit('Exiting in dftg_cdfit_gen.',-9)
      do 305 ishp = ishatf, ishatl
c!!!      do 305 ishp = 1, nshells_cd
c
         if (ishp.eq.next)then
c
            next = nxtask(nproc,igran) + 1
c
            if (.not. bas_cn2bfr(CD_bas_han,ishp,ifirstp,ilastp))
     &         call errquit('Exiting in dftg_cdfit_gen.',9)
               nshbfp = ilastp - ifirstp + 1
c
c!!!            do ishq = 1, ishp
            do ishq = 1, ishatf-1
c
               if (.not. bas_cn2bfr(CD_bas_han,ishq,ifirstq,ilastq))
     &            call errquit('Exiting in dftg_cdfit_gen.',10)
               nshbfq = ilastq - ifirstq + 1
               Nintegrals = nshbfp*nshbfq
c
               if (do_force.or.do_cpks_r) then
                  call intd_2e2c(CD_bas_han, ishp, CD_bas_han, ishq, 
     &                 lscr, scr, lbuf, buf, idatom)
               else if (do_hess) then
                  call intdd_2e2c(CD_bas_han, ishp, CD_bas_han, ishq, 
     &                 lscr, scr, lbuf, buf, idatom)
               endif
c
               fac = -1.d0
               if (ishp.eq.ishq)fac=-.5d0
c
               if (do_force) then
c
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
c
               else if (do_hess) then
c
                  write(*,*)'Are dimensions right on p and q?'
                  call d2int_2e2c_block_jhess(CD_coef(ifirstp),
     &                 CD_coef(ifirstq), hess, nat, buf, fac, idatom,
     &                 nshbfp, nshbfq)
c
               else if (do_cpks_r) then
c
               write(*,*)'Are dimensions right on p and q?'
                  call dint_2e2c_block_jdmat(CD_coef, db, nbf_cd, buf,
     &                 idatom, ifirstp, ilastp, ifirstq, ilastq)
c
               endif
            enddo
         endif
 305  continue
 304  continue
c
      next=nxtask(-nproc,igran)
c
      if (do_cpks_r) then
         write(*,*)'*** Do global sum of db'
      endif
c
      return 
      end 
