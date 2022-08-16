      Subroutine dftg_cdfit(geom,AO_bas_han, CD_bas_han,
     &                      nbf_cd, nat, tol2e, scr, 
     &                      lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                      iga_dens, force, CD_coef,
     &                      oskel, nfit, adft)
c
c     Turned into wrapper to call generalized function
c
c     BGJ - 3/99
c    
*
* $Id$
*
      integer geom
      integer nfit
      integer  AO_bas_han, CD_bas_han
      integer nbf_cd, nat, max_sh_bf
      integer iga_dens(*), lbuf, lscr
      double precision scr(lscr), buf(lbuf)
      double precision vec(*),
     &                 PPP(*)
!     &                 PPP(max_sh_bf,max_sh_bf)
      double precision CD_coef(nbf_cd,3*nfit)
      double precision force(3,nat)
      double precision tol2e
      logical oskel, adft, doit
c
      call dftg_cdfit_gen(geom, AO_bas_han, CD_bas_han,
     &                    nbf_cd, nat, tol2e, scr, 
     &                    lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                    iga_dens, force, 0d0, 0, 0d0,
     &                    CD_coef, oskel, 1, nfit, adft)
c
      return
      end
c
c     Generalized from dftg_cdfit to calculate fitted J derivative
c     matrices and hessian as well as fitted J forces
c
c     BGJ - 3/99
c
      Subroutine dftg_cdfit_gen(geom, AO_bas_han, CD_bas_han,
     &                           nbf_cd, nat, tol2e, scr, 
     &                          lscr, buf, lbuf, vec, PPP, max_sh_bf,
     &                          iga_dens, force, hess,
     &                          g_dj, db, CD_coef, oskel_in, calc_type,
     &                          nfit, adft)
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "sym.fh"
#include "eaf.fh"
#include "stdio.fh"
#include "inp.fh"
c
      integer nfit
      integer geom
      integer  AO_bas_han, CD_bas_han
      integer  nbf_cd, nat, max_sh_bf
      integer iga_dens(*), lbuf, lscr, calc_type
c
c     scratch arrays for calls to integrals package
c
      double precision scr(lscr), buf(lbuf)
      integer idatom(4)
c
c     scratch arrays for dgem
c
      double precision vec(*),
     &                 PPP(*)
c
c-->  Charge Density Expansion Coefficients
c     
      double precision CD_coef(nbf_cd,3*nfit)
                       ! nfit = 1:
                       ! cd_coef(:,1) -> CD_tot
                       ! cd_coef(:,2) -> XC_tot
                       ! cd_coef(:,3) -> Work_tot

                       ! nfit = 2:
                       ! cd_coef(:,1) -> CD_up
                       ! cd_coef(:,2) -> CD_dw
                       ! cd_coef(:,3) -> XC_up
                       ! cd_coef(:,4) -> XC_dw
                       ! cd_coef(:,5) -> Work_up
                       ! cd_coef(:,6) -> Work_dw
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
      logical doesit, oskel, adft
      integer nproc,me,icount,igran
      integer nshells_ao,nshells_cd,ishp,ishq
      integer ifirstq, ilastq, ifirstp, ilastp,nshbfq,nshbfp,Nao2
      integer Nintegrals,ishd,ifirstd,ilastd,nshbfd
      integer iat,next,nxtask,icart,l_dfkl,k_dfkl
      integer iatom,ishatf,ishatl
      integer atom_p,atom_q
      integer sh_lo_p,sh_hi_p, sh_lo_q, sh_hi_q
      integer bf_lo_p,bf_hi_p,bf_lo_q,bf_hi_q
      integer k_p,l_p
      integer nbf_at_p
      integer icount_p,icount_q,nbf_at_q
      integer ishq_0,ishq_1,icount_task
      integer nbf_ao_mxnbf_ce
      double precision PPP_max,dE,fac,schwarz_shell,q1
      double precision offset,cd_max,cd_maxshell
      logical do_force, do_cpks_r, do_hess, doit
      integer ierr,LU
      integer k_at,l_at,atom_p_in,atom_q_in,nat_used
      integer ifit,lencd
      character*255 errmsg
c
      character*(nw_max_path_len) cdfit
c
      double precision time0
      double precision gradient
      double precision dabsmax
      double precision screen
      external nxtask,schwarz_shell,dabsmax
      double precision ddot
      external ddot
      nproc  = ga_nnodes()
      me = ga_nodeid()
c
      do_force = calc_type .eq. 1
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_force .or. do_cpks_r .or. do_hess)) then
         call errquit('dftg_cdfit_gen: illegal calculation type',1,
     &       INPUT_ERR)
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
     &        errquit('dftg_cdfit_gen: cannot allocate dfkl',1,
     &       MA_ERR)
         call dfill(nbf_cd*3*nat, 0.0d0, db, 1)
      endif
c
c     read CD_coeff vector written by fitcd
c
      if (me.eq.0)then
         call util_file_name('cdfit', .true., .false., cdfit)
         inquire(file=cdfit,exist=doesit)
         if (.not. doesit)call errquit
     &        ('dftg_cdfit_gen: file cdfit does not exist',0, DISK_ERR)
         LU=79
         ierr=eaf_open(cdfit(1:inp_strlen(cdfit)),eaf_rw, LU)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('dftg_fitcd: eaf_open failed', ierr, DISK_ERR)
        endif
        do ifit=1,nfit
          offset = (ifit-1)*nbf_cd*8d0
          ierr = eaf_read(LU, offset, cd_coef(1,ifit), 8*nbf_cd)
          if (ierr .ne. 0) then
             call eaf_errmsg(ierr, errmsg)
             write(LuOut,*) ' IO error message ',
     $            errmsg(1:inp_strlen(errmsg))
             call errquit('dftg_fitcd: read failed',0, DISK_ERR)
          endif
        enddo

        ! Read XC fitting coefficients as well
        if (adft) then
          do ifit=nfit+1,2*nfit
            offset = (ifit-1)*nbf_cd*8d0
            ierr = eaf_read(LU, offset, cd_coef(1,ifit), 8*nbf_cd)
            if (ierr .ne. 0) then
               call eaf_errmsg(ierr, errmsg)
               write(LuOut,*) ' IO error message ',
     $              errmsg(1:inp_strlen(errmsg))
               call errquit('dftg_fitcd: read failed',0, DISK_ERR)
            endif
          enddo
        endif
        ierr=eaf_close(LU)
        if (ierr .ne. 0) then
           write(luout,*) ga_nodeid(),' closing FD =',LU,
     ,          ' cdfit '
           call eaf_errmsg(ierr, errmsg)
           write(LuOut,*) ' IO error message ',
     $          errmsg(1:inp_strlen(errmsg))
           call util_flush(luout)
           call errquit('dftg_fitcd: eaf_close failed',0, DISK_ERR)
        endif
      endif
c
      lencd = nbf_cd*nfit
      if (adft) lencd = 2*lencd

      call ga_brdcst(Msg_dftg_cdf,CD_coef,
     M     ma_sizeof(MT_DBL,lencd,MT_BYTE),0)
c
      ! Get CD_tot in working entry 2*nfit+1
      call dcopy(nbf_cd,cd_coef(1,1),1,cd_coef(1,2*nfit+1),1)
      if (nfit.gt.1) then
        call daxpy(nbf_cd,1d0,cd_coef(1,2),1,cd_coef(1,2*nfit+1),1)
      endif

      ! For ADFT we need to XC coeff as well
      if (adft) then
        ! Open shell cases must be spin-resolved, so first copy
        ! CD_tot to last entry (3*nfit) and then add XC_dw
        if (nfit.gt.1) then
          call dcopy(nbf_cd,cd_coef(1,2*nfit+1),1,cd_coef(1,3*nfit),1)
          call daxpy(nbf_cd,1d0,cd_coef(1,2*nfit),1,cd_coef(1,3*nfit),1)
        endif
        ! Add XC_up to entry 2*nfit+1
        call daxpy(nbf_cd,1d0,cd_coef(1,nfit+1),1,cd_coef(1,2*nfit+1),1)
      endif

      cd_max = dabsmax(nbf_cd,cd_coef(1,2*nfit+1))
      if (adft) then
        cd_max = max(cd_max,dabsmax(nbf_cd,cd_coef(1,3*nfit)))
      endif
c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
      if (.not. bas_numcont(AO_bas_han,nshells_ao))then
         call errquit('Exiting in dftg_cdfit_gen.',2, BASIS_ERR)
      endif
c      
      if (.not. bas_numcont(CD_bas_han,nshells_cd))then
         call errquit('Exiting in dftg_cdfit_gen.',3, BASIS_ERR)
      endif
c
c     2el. 3-c integral derivatives
c
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     C     call errquit('dftg_fitcd: bas_nbf failed ',0,BASIS_ERR)
      if(.not.MA_Push_Get(MT_Dbl,nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce*nfit,
     G     'gather vector',l_p,k_p))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
      if(.not.MA_Push_Get(MT_int,nat,'at list',l_at,k_at))
     &  call errquit('dft_cdfit: cannot allocate at list',0, MA_ERR)
c     sort atoms list based on decreasing atomic weight
      nat_used=nat
      call dft_sortatbyz_screen(ao_bas_han,geom,nat_used,int_mb(k_at))
c
      next = nxtask(nproc,1) + 1
      icount_task=1
c
      do atom_p_in=1,nat_used 
         atom_p=int_mb(k_at+atom_p_in-1)
        if (.not. bas_ce2cnr( AO_bas_han, atom_p, sh_lo_p, sh_hi_p)) 
     &       call errquit('Exiting in dftg_cdfit',110, BASIS_ERR)
        if (.not. bas_ce2bfr(AO_bas_han,atom_p,bf_lo_p,bf_hi_p)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        nbf_at_p = bf_hi_p - bf_lo_p + 1
        do atom_q_in=atom_p_in, nat_used
         if (icount_task.eq.next)then
            atom_q=int_mb(k_at+atom_q_in-1)
           if (.not. bas_ce2cnr(AO_bas_han,atom_q,sh_lo_q,sh_hi_q)) 
     &       call errquit('Exiting in dftg_cdfit',110, BASIS_ERR)
           if (.not. bas_ce2bfr(AO_bas_han,atom_q,bf_lo_q,bf_hi_q)) 
     &          call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
           nbf_at_q=bf_hi_q - bf_lo_q + 1

               call ga_get(iga_dens(1), bf_lo_q, bf_hi_q, 
     C                bf_lo_p,bf_hi_p,dbl_mb(k_p),nbf_at_q)

               if (adft .and. nfit.gt.1) then
                 call ga_get(iga_dens(2), bf_lo_q, bf_hi_q, 
     C                bf_lo_p,bf_hi_p,
     C                dbl_mb(k_p+nbf_at_q*nbf_at_p),nbf_at_q)
               endif

               icount_p=0
           do  ishp = sh_lo_p, sh_hi_p
c
c
            if (.not. bas_cn2bfr(AO_bas_han,ishp,ifirstp,ilastp))
     &         call errquit('Exiting in dftg_cdfit_gen.',6, BASIS_ERR)
            nshbfp = ilastp - ifirstp + 1
c     
               ishq_0 =  sh_lo_q
               ishq_1 =  sh_hi_q
cicic               if(atom_p.eq.atom_q) ishq_0=ishp
               icount_q=0
            do ishq = ishq_0,ishq_1
               if (.not. bas_cn2bfr(AO_bas_han,ishq,ifirstq,ilastq))
     &            call errquit('Exiting in dftg_cdfit_gen.',7,
     &       BASIS_ERR)
               nshbfq = ilastq - ifirstq + 1
               Nao2 = nshbfp*nshbfq
c
               if (do_force.or.do_hess) then
                  fac = 2d0
c                  if (ishp.eq.ishq) fac = 1d0
                  if (atom_p.eq.atom_q) fac = 1d0
               else
                  fac = 1d0
               endif
c
c              get shell block of DM
c

c     copy pmat into ppp for now
                  call dft_copyfrom(
     I                   nbf_at_q,nshbfq,nshbfp,
     I                   dbl_mb(k_p),ppp,icount_q,icount_p)

                  if (adft .and. nfit.gt.1) then
                    call dft_copyfrom(
     I                   nbf_at_q,nshbfq,nshbfp,
     I                   dbl_mb(k_p+nbf_at_p*nbf_at_q),ppp(nao2+1),
     .                   icount_q,icount_p)
                  endif
               
c
c              Schwarz screening on product MAX(ABS(DM(i,j)))*schwarz(i,j)
c
               PPP_max = dabsmax(nao2,PPP)
               if (adft .and. nfit.gt.1) then
                 PPP_max = max(PPP_max,dabsmax(nao2,PPP(nao2+1)))
               endif

               screen = tol2e/PPP_max/schwarz_shell(ishp,ishq)
               if (cd_max.gt.screen)then
                 
                  do ishd = 1, nshells_cd
                     if (.not. bas_cn2bfr(CD_bas_han, ishd, 
     &                  ifirstd, ilastd))
     &                  call errquit('Exiting in dftg_cdfit_gen.',8,
     &       BASIS_ERR)
                     nshbfd = ilastd - ifirstd + 1
                     Nintegrals = nshbfp*nshbfq*nshbfd
                     doit = dabsmax(nshbfd,cd_coef(ifirstd,2*nfit+1))
     &                      .gt.screen
                     if (adft.and.(nfit.gt.1)) then
                       doit = doit .and.
     &                 dabsmax(nshbfd,cd_coef(ifirstd,3*nfit)).gt.screen
                     endif

                     if (.not.doit) cycle
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
                           call dgemv('t',nao2,3*nshbfd,1d0,
     &                       buf(icount),nao2,ppp,1,0d0,vec,1)

                           ! For standard CDFIT or (ADFT.and.(ipol.eq.1))
                           ! this step is
                           !     Ptot * d3cERI * (CD_tot + XC_tot)
                           !
                           ! for ADFT and ipol > 1:
                           !     P_up * d3cERI * (CD_tot + XC_up)

                           do icart=1,3
                             dE = dot_product(vec((icart-1)*nshbfd+1:
     &                                        icart*nshbfd),
     &                            cd_coef(ifirstd:ilastd,2*nfit+1))
                             force(icart,idatom(iat)) = 
     &                       force(icart,idatom(iat)) + dE * fac * q1
                           enddo

                           ! Need to complete ADFT open-shell case with
                           !    P_dw * d3cERI * (CD_tot + XC_dw)
                           !
                           ! so we end up with
                           !    P_tot * d3cERI * CD_tot +
                           !    P_up * d3cERI * XC_up +
                           !    P_dw * d3cERI * XC_dw
                           if (adft.and.(nfit.eq.2)) then
                             call dgemv('t',nao2,3*nshbfd,1d0,
     &                         buf(icount),nao2,ppp(nao2+1),1,0d0,vec,1)
                             do icart=1,3
                               dE = dot_product(vec((icart-1)*nshbfd+1:
     &                                        icart*nshbfd),
     &                            cd_coef(ifirstd:ilastd,3*nfit))
                               force(icart,idatom(iat)) = 
     &                         force(icart,idatom(iat)) + dE * fac * q1
                             enddo
                           endif
                         endif
                         icount = icount + 3*Nintegrals
                     enddo
c
                     else if (do_hess) then
c
                        write(*,*)'Are dimensions right on ints and P?'
                        call d2int_2e3c_block_jhess(
     &                       CD_coef(ifirstd,2*nfit+1),
     &                       PPP, nshbfq, hess, nat, buf, fac*q1,
     &                       idatom, nshbfd, nshbfp, nshbfq)
c
                     else if (do_cpks_r) then
c
                        write(*,*)'Are dimensions right on ints and P?'
                        call dint_2e3c_block_jdmat(
     &                       CD_coef(ifirstd,2*nfit+1),
     &                       PPP, nbf_cd, nshbfq, g_dj,
     &                       db(ifirstd,1,1), dbl_mb(k_dfkl), buf,
     &                       ishp .eq. ishq, idatom, nshbfd,
     &                       ifirstp, ilastp, ifirstq, ilastq)
c
                     endif
  204                continue
                  enddo
               endif
               icount_q=icount_q+nshbfq
            enddo ! ishq
            icount_p=icount_p+nshbfp
      enddo                     ! ishp
      next = nxtask(nproc,1) + 1
      endif
      icount_task=icount_task+1
      enddo                     ! atom_q
      enddo                     ! atom_p
c
      next = nxtask(-nproc,1)
c
      if (.not.MA_chop_stack(l_p))
     &     call errquit('dftg_cdfit_gen: cannot pop l_p',1, MA_ERR)
      if (do_cpks_r) then
         if (.not.MA_pop_stack(l_dfkl))
     &        call errquit('dftg_cdfit_gen: cannot pop dfkl',1, MA_ERR)
      endif
      if (do_cpks_r) then
         write(*,*)'DB after 3c part:'
         call output(db, 1, nbf_cd, 1, 3*nat, nbf_cd, 3*nat, 1)
c!!!         call dfill(nbf_cd*3*nat, 0.0d0, db, 1)
      endif


      if(.not.MA_Push_Get(MT_Dbl,lbuf,'gather vector',l_p,k_p))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
c
c     2el. 2-c integral derivatives
c
      igran = (nshells_ao + mod(nshells_ao,2))/(4*nproc)
      igran = max(1,igran)
c
      next = nxtask(nproc,igran) + 1
c
c     Since 2e2c integral derivatives are zero unless the centers
c     involved are different, the loops are set up to run only over
c     pairs of shells on different centers, saving a small amount
c     of work
c
      ! Prepare the fitting coefficients for the 2cERI part in AFDT.
      ! Here, we need our working fields to be
      !
      ! Closed-shell
      !   cd_coef(:,3) -> CD_tot + 2*XC_tot
      !
      ! Open-shell
      !   cd_coef(:,5) -> CD_tot + 2*XC_up
      !   cd_coef(:,6) -> CD_tot + 2*XC_dw
      !
      ! (see bldcblk below)
      !
      ! We already the XC_coeffs included in the working fields, 
      ! so we only need to add again the corresponding XC coeffs.
      if (adft) then
        call daxpy(nbf_cd,1d0,cd_coef(1,nfit+1),1,cd_coef(1,2*nfit+1),1)
        if (nfit.gt.1) then
          call daxpy(nbf_cd,1d0,cd_coef(1,2*nfit),1,cd_coef(1,3*nfit),1)
        endif

      ! We will use cd_coef(:,1) to handle CD_tot in standard 
      ! CD fitting
      elseif(nfit.gt.1)then
        call daxpy(nbf_cd,1d0,cd_coef(1,2),1,cd_coef(1,1),1)
      endif

      do 304 iatom = 1, nat
         if (.not. bas_ce2cnr(CD_bas_han,iatom,ishatf,ishatl))
     &        call errquit('Exiting in dftg_cdfit_gen.',-9, BASIS_ERR)
      do 305 ishp = ishatf, ishatl
c
         if (ishp.eq.next)then
c
            next = nxtask(nproc,igran) + 1
c
            if (.not. bas_cn2bfr(CD_bas_han,ishp,ifirstp,ilastp))
     &         call errquit('Exiting in dftg_cdfit_gen.',9, BASIS_ERR)
               nshbfp = ilastp - ifirstp + 1
c
            do ishq = 1, ishatf-1
c
               if (.not. bas_cn2bfr(CD_bas_han,ishq,ifirstq,ilastq))
     &            call errquit('Exiting in dftg_cdfit_gen.',10,
     &       BASIS_ERR)
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
               call dcopy(nshbfq*nshbfp,0d0,0,dbl_mb(k_p),1)
               call bldcblk(nshbfq,nshbfp,ifirstq,ifirstp,
     $                      dbl_mb(k_p),cd_coef,cd_coef(1,2*nfit+1))
               if (adft) then
                 if (nfit.gt.1) then
                   call bldcblk(nshbfq,nshbfp,ifirstq,ifirstp,
     $                      dbl_mb(k_p),cd_coef(1,2),cd_coef(1,3*nfit))
                 endif
               endif

               icount = 1
               do iat = 1, 2
                  if (idatom(iat).gt.0)then
                    call dgemv('t',Nintegrals,3,fac,buf(icount),
     $                    Nintegrals,dbl_mb(k_p),1,1.0d0,
     $                    force(1,idatom(iat)),1)
                  endif
                  icount = icount + 3*Nintegrals
               enddo
c
               else if (do_hess) then
c
                  write(*,*)'Are dimensions right on p and q?'
                  call d2int_2e2c_block_jhess(CD_coef(ifirstp,1),
     &                 CD_coef(ifirstq,1), hess, nat, buf, fac, idatom,
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

      if (.not.MA_chop_stack(l_p))
     &     call errquit('dftg_cdfit_gen: cannot pop l_p',1, MA_ERR)
c
      return 
      end 

      subroutine bldcblk(na,nb,inia,inib,buf,cd,cdxc)
      ! We need to contract the derivatives of the 2cERIs
      ! we the appropriate fitting coefficients.
      !
      ! This amounts to build to block: 
      !
      ! blk = 0.5*CD_ispin(i)*CD_tot(j) + CD_ispin(i)*XC_ispin(j)
      !     = 0.5*CD_ispin(i)*(CD_tot(j) + 2*XC_ispin(j))
      !     = 0.5*CD_ispin(i)*CDXC_ispin(j)
      !
      ! (The factor 2 was included earlier)
      !
      ! For closed-shells and standard CDFIT, this subroutine will
      ! be called only once with CD_tot in the CD_ispin fields,
      ! and XC_tot in the XC_ispin field.
      ! 
      ! For open-shell ADFT, this subroutine will be called twice
      ! to get
      !
      ! 0.5*CD_tot*d2cERI*CD_tot + CD_up*d2cERI*XC_up + CD_dw*d2cERI*XC_dw
      !

      implicit none
      integer na,nb,inia,inib
      double precision buf(na,nb),cd(*),cdxc(*)

      integer i,j
      do j=1,nb
        do i=1,na
          buf(i,j) = buf(i,j) + 0.5d0*(cd(i+inia-1)*cdxc(j+inib-1) +
     $               cdxc(i+inia-1)*cd(j+inib-1))
        enddo
      enddo
      end
