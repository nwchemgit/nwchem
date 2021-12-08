      subroutine tddft_fock_fit(nder,ndens,jfac,tol2e,g_dens,g_fock) 

      implicit none
#include "errquit.fh"
c
#include "cdft.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "schwarz.fh"
#include "sym.fh"

      integer ndens
      integer nder
      integer g_dens
      integer g_fock
      integer tol2e
      double precision jfac(*)

      integer lo(3), hi(3), ld(2)
      integer ilo, ihi, jlo, jhi, adrc, ldc
      integer nproc, me, n1, i1

      logical dowork, v_nonzero, doshell, iolgc

      integer nao2
      integer lTvec, iTvec
      integer l_p, k_p
      integer lgvec, igvec
      integer leri, ieri
      integer lscr, iscr, nscr
      integer lcdc, icdc
      integer lcd_shell, icd_shell
      integer maxg, mscratch_2e3c

      integer atom_c,ishc,sh_lo_c,sh_hi_c,nbf_at_c,nshbfc,nshbf_cd
      integer bf_lo_c,bf_hi_c
      integer bf_lo_d,bf_hi_d
      integer atom_d,ishd,sh_lo_d,sh_hi_d,nbf_at_d,nshbfd
      integer ishd_0, ishd_1, ifirstc, ifirstd, ilastc, ilastd
      integer atom_p,ishp,ifirstp,ilastp,nshbfp
      integer size_atblock
      integer lmiss,iptr,j
      integer ishbf,lshbf

      integer icount, icount_c, icount_d
      integer istart, nintegrals

      integer idens,i,nnii,nfunz,ninteg_api

      integer g_tmp

      double precision threshold, az, shmax, cd_coef_max
      double precision pmat_max,eri_est,fac,vcoul_est,q1,qscale,qsym

      integer,external :: nxtask_static, nxtask
      double precision, external :: ddot, dabsmax

      me = ga_nodeid()
      iolgc = noio.ne.1

      nproc  = ga_nnodes()
      THRESHOLD = 10.d0**(-itol2e)

      call ga_zero(g_fock)

      if (sum(abs(jfac(1:ndens))).eq.0d0) goto 100


      if(.not.MA_Push_Get(MT_Dbl,nbf_cd*ndens,'Tvec',lTvec,iTvec))
     &  call errquit('dft_fitcd: cannot allocate Tvec',0, MA_ERR)
      call dcopy(nbf_cd*ndens, 0.D0, 0, Dbl_MB(iTvec), 1)
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd*ndens,'CDCvec',lcdc,icdc))
     &  call errquit('dft_fitcd: cannot allocate CDCvec',0, MA_ERR)
      call dcopy(nbf_cd*ndens, 0.D0, 0, Dbl_MB(icdc), 1)

      Nao2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(.not.MA_Push_Get(MT_Dbl,Nao2*ndens,'gather vector',l_p,k_p))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
      Nao2 = nbf_ao_mxnbf_cn*nbf_ao_mxnbf_cn
      if(.not.MA_Push_Get(MT_Dbl,Nao2*ndens,'gathervector',
     &        lgvec,igvec))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)

      call int_mem_2e3c(maxg, mscratch_2e3c)
      if(.not.MA_Push_Get(MT_Dbl,maxg*nshells_cd,'ERI',lERI,iERI))
     &  call errquit('dft_fitcd: cannot allocate ERI',0, MA_ERR)
      nscr =  mscratch_2e3c 
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0, MA_ERR)

      icount = 0
      n1 = nxtask_static(nproc,1)

      do atom_c=1,ncenters
        if (.not. bas_ce2cnr( AO_bas_han, atom_c, sh_lo_c, sh_hi_c)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        if (.not. bas_ce2bfr( AO_bas_han, atom_c, bf_lo_c, bf_hi_c)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        nbf_at_c=bf_hi_c - bf_lo_c + 1

        do atom_d=1,atom_c
          dowork=n1.eq.icount
          if(dowork) then

            if (.not.bas_ce2bfr(AO_bas_han,atom_d,bf_lo_d,bf_hi_d)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
            if (.not.bas_ce2cnr(AO_bas_han,atom_d,sh_lo_d,sh_hi_d)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
            nbf_at_d=bf_hi_d - bf_lo_d + 1

c     get DM atom block
            size_atblock= nbf_at_d*nbf_at_c
            do idens = 1, ndens
               lo(1:3) = (/ idens, bf_lo_d, bf_lo_c /)
               hi(1:3) = (/ idens, bf_hi_d, bf_hi_c /)
               ld(1:2) = (/ 1, bf_hi_d - bf_lo_d + 1 /)
               call nga_get(g_dens,lo,hi,
     .              dbl_mb(k_p+(idens-1)*size_atblock),ld)
            enddo
            icount_c=0
            do ishc = sh_lo_c,sh_hi_c
               if (.not. bas_cn2bfr(AO_bas_han,ishc,ifirstc,ilastc)) 
     &          call errquit('Exiting in dft_fitcd',111, BASIS_ERR)
               nshbfc = ilastc-ifirstc+1

               icount_d=0
               do ishd =  sh_lo_d, sh_hi_d
c
                 if (.not. bas_cn2bfr(AO_bas_han,ishd,ifirstd,ilastd))
     &             call errquit('Exiting in dft_fitcd',1111, BASIS_ERR)
                 nshbfd = ilastd-ifirstd+1
c
c                Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c
                 ERI_est = schwarz_shell(ishc,ishd)
                 if(ERI_est .gt.THRESHOLD) then
                   Nao2 = nshbfd*nshbfc
                   pmat_max=1d0

                   do idens = 1, ndens
                     call dft_copyfrom(
     I                     nbf_at_d,nshbfd,nshbfc,
     I                     dbl_mb(k_p+(idens-1)*size_atblock),
     $                     dbl_mb(igvec+(idens-1)*Nao2),
     A                     icount_d,icount_c)
                   enddo
c     
c                 Determine the maximum density matrix element for this product of AO shells.
c     
                   Pmat_max = dabsmax(Nao2*ndens,DBL_MB(igvec))
                   if( Pmat_max*ERI_est.gt.THRESHOLD*1d-3 )then

                     istart = iERI
                     call dcopy(nao2*nbf_cd,0d0,0,dbl_mb(istart),1)

                     fac = 2.0d0
                     if(atom_c.eq.atom_d) fac = 1.d0

                     do ishp = 1,nshells_cd
c                  
c                   Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c                  
                       if (.not.bas_cn2bfr(CD_bas_han,ishp,ifirstp,
     .                                    ilastp))
     &                    call errquit('fitcd.',5, BASIS_ERR)
                       nshbfp = ilastp-ifirstp+1
c                  
                       nintegrals = nshbfp*Nao2
                       Ninteg_api = maxg
c
                       q1=1.d0
                       doshell=.true.
                       if (doshell) then
                          call int_2e3c(CD_bas_han,ishp,AO_bas_han,
     I                             ishc,ishd,nscr,Dbl_MB(iscr),
     N                             Ninteg_api,Dbl_MB(istart))
                          call dscal(nintegrals,q1,dbl_mb(istart),1)
                       else
                          call dcopy(nintegrals,0d0,0,dbl_mb(istart),1)
                       endif
                       istart = istart+nintegrals
c     
                     enddo       !ishp
                     do idens=1,ndens
                       call dgemv('T',Nao2,nbf_cd,fac,dbl_mb(ieri),
     .                             nao2,dbl_mb(igvec+(idens-1)*nao2),1,
     T                            1.d0,dbl_mb(iTvec+(idens-1)*nbf_cd),1)
                     enddo
                   endif
                 endif
                 icount_d = icount_d + nshbfd
              enddo
              icount_c = icount_c  + nshbfc
            enddo
            n1 = nxtask_static(nproc,1)
          endif
          icount = icount + 1
        enddo
      enddo

      n1 = nxtask_static(-nproc,1)
c
      call ga_mask_sync(.false.,.true.)
      call ga_dgop(1,dbl_mb(iTvec),nbf_cd*ndens,'+')
c
      call ga_sync()

      if (IOLGC)then
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERInv', 
     &                       nbf_cd, 0, g_cdinv)) 
     &      call errquit('dft_fitcd: error creating g_cdinv',0, GA_ERR)
        if (me.eq.0)then
           lmiss = 0
           call dft_invio('CDI',g_cdinv,nbf_cd,'READ',lmiss)
        endif
        call ga_brdcst(2000, lmiss,MA_sizeof(MT_INT,1,MT_BYTE),0)
        if(lmiss.ne.0)call errquit
     &       (' dft_fitcd; dft_invio - abnormal read of CDI ',0, GA_ERR)
      endif
c
c     Multiply T and cgtf by S**(-1) and store
c
      call  ga_distribution(g_cdinv,me,ilo,ihi,jlo,jhi)
      if (ilo.gt.0) then
        call ga_access(g_cdinv, ilo, ihi, jlo, jhi, adrc, ldc)
        nnii=ihi-ilo+1
        do j = jlo,jhi
          iptr=adrc+(j-jlo)*ldc
          do idens = 1, ndens
             AZ = ddot(nnii,DBL_mb(iptr),1,
     &            DBL_MB(iTvec+(idens-1)*nbf_cd+ilo-1),1)
             dbl_mb(icdc+(idens-1)*nbf_cd + j - 1) = az        
          enddo
        enddo
        call ga_release(g_cdinv, ilo, ihi, jlo, jhi)
      endif! own patch
      call ga_mask_sync(.true.,.false.)
      call ga_dgop(1,dbl_mb(icdc),nbf_cd*ndens,'+')
      if (.not.ga_destroy(g_cdinv))
     $  call errquit("tddft_fock_fit: destroy GA",g_cdinv,GA_ERR) 

      if(.not.MA_Push_Get(MT_Dbl,nshells_cd,'CD shell max',
     &  lcd_shell,icd_shell))
     &  call errquit('dft_fitvc: cannot allocate CD shell max',0,0)
c
      CD_coef_max = 0.D0
      i1 = 1
      do  ishp = 1,nshells_cd
        if( .not. bas_cn2bfr(CD_bas_han,ishp,ishbf,lshbf) )then
          call errquit('Exiting in fitvcoul.',4,0)
        end if
        shmax = 0.D0
        do i = ishbf,lshbf
           do idens = 1, ndens
              shmax = max(abs(dbl_mb(icdc+(idens-1)*nbf_cd+i-1)), shmax)
           enddo
        enddo
        Dbl_MB((icd_shell-1)+ishp) = shmax
        CD_coef_max = max(shmax, CD_coef_max)
      enddo

      icount = 0
      n1 = nxtask(nproc,1)
      do atom_c = 1,ncenters
         if (.not. bas_ce2cnr( AO_bas_han, atom_c, sh_lo_c, sh_hi_c)) 
     &        call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
         if (.not. bas_ce2bfr( AO_bas_han, atom_c, bf_lo_c, bf_hi_c)) 
     &        call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
         nbf_at_c=bf_hi_c - bf_lo_c + 1
c     
         do atom_d= 1,atom_c
            dowork=n1.eq.icount
            if(dowork) then
               v_nonzero=.false.
               if (.not.bas_ce2bfr(AO_bas_han,atom_d,bf_lo_d,bf_hi_d)) 
     &              call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
               nbf_at_d=bf_hi_d - bf_lo_d + 1
               if (.not.bas_ce2cnr(AO_bas_han,atom_d,sh_lo_d,sh_hi_d))
     &              call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
               icount_c=0
czero k_p
               call dcopy(nbf_at_d*nbf_at_c*ndens,0d0,0,dbl_mb(k_p),1)
               do ishc=sh_lo_c,sh_hi_c

                  if( .not. bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &                 call errquit('Exiting in fitvcoul.',3,0)
                  nshbfc=ilastc-ifirstc+1
                  
                  ishd_0 = sh_lo_d
                  ishd_1 = sh_hi_d
                  if (atom_c.eq.atom_d) ishd_1=ishc

                  icount_d=0
                  do ishd =  ishd_0,  ishd_1
                     
                     if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, 
     I                    ilastd))
     &                    call errquit('Exiting in fitvcoul.',4,0)
                     nshbfd=ilastd-ifirstd+1
                     Nao2 = nshbfc*nshbfd
c           
c     Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c           
                     ERI_est=schwarz_shell(ishc,ishd)
                     if( ERI_est*CD_coef_max.gt.THRESHOLD )then

                        v_nonzero=.true.
                        qscale=1d0
                        if(ishc.ne.ishd) qscale=2d0
                           call dcopy(Nao2*ndens,0.D0,0,dbl_mb(igvec),1)
                           nfunz=0
                           do ishp = 1,nshells_cd
c     
                              if( .not.bas_cn2bfr(CD_bas_han, ishp, 
     I                             ifirstp, ilastp))
     &                             call errquit('fitvc',0,BASIS_ERR)
                              nshbf_cd=ilastp-ifirstp+1
                              Vcoul_est=dbl_mb(icd_shell+ishp-1)*ERI_est
                              if( Vcoul_est.gt.THRESHOLD )then
                                 Ninteg_api=maxg
c                  
c     Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c     
c     check for use of symmetry
                                 qsym=1d0
                                 doshell=.true.

                                 if(doshell) then
                                    call int_2e3c(CD_bas_han,ishp,
     A                                   AO_bas_han,ishc,ishd,nscr,
     D                                   Dbl_MB(iscr),Ninteg_api,
     &                                   Dbl_MB(iERI))
c     
c     Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c     
                                    q1=qscale*qsym
                                    do idens=1,ndens
                                    call dgemv('N',Nao2,nshbf_cd,q1,
     &                                   Dbl_MB(iERI),Nao2,
     .                             dbl_mb(icdc+nfunz+(idens-1)*nbf_cd),
     &                                   1,1.D0,
     &                                   Dbl_MB(igvec+(idens-1)*nao2),1)
                                    enddo
                                 endif
                              endif ! schw on VCoul_est
c
                              nfunz=nfunz+nshbf_cd
c               
                           enddo ! ishp
c
c     Scatter the computed values of Vcoul.
c           

                           do idens=1,ndens
                           call dft_copyto(
     I                          nbf_at_d,nshbfd,nshbfc,
     I                          dbl_mb(k_p+(idens-1)*nbf_at_c*nbf_at_d),
     .                          dbl_mb(igvec+(idens-1)*nao2),
     A                          icount_d,icount_c)
                           enddo
                        endif   ! schwartz with cdmax

                     icount_d=icount_d+nshbfd
                  enddo         ! ishd
                  icount_c=icount_c+nshbfc
               enddo            ! ishc
c     ga_put only if anything ever survived the schwarz scr
               
               if(v_nonzero) then
                 do idens=1,ndens
                    lo(1:3) = (/ idens, bf_lo_d, bf_lo_c /)
                    hi(1:3) = (/ idens, bf_hi_d, bf_hi_c /)
                    ld(1:2) = (/ 1, nbf_at_d /)
                    call nga_put(g_fock, lo, hi,
     B                dbl_mb(k_p+(idens-1)*nbf_at_d*nbf_at_c),ld)
                 enddo
               endif

               n1 = nxtask(nproc,1)
            endif
c     
          icount = icount + 1
c     
       enddo                    ! atom_d
      enddo                     ! atom_c
c
c      lo(2) = 1
c      lo(3) = 1
c      hi(2) = nbf
c      hi(3) = nbf
c      do idens=1,ndens
c         lo(1)=idens
c         hi(1)=idens
c         call nga_print_patch(g_fock,lo,hi,0)
c      enddo
c
      if(.not.ma_chop_stack(ltvec))
     &  call errquit('dft_fitcd: cannot pop stack',12, MA_ERR)
      n1 = nxtask(-nproc,1)

  100 continue

      call xc_newfock(geom,ao_bas_han,jfac,jfac,nbf_ao,
     &        nder, ndens, ndens, g_fock,g_dens)

      call ga_mask_sync(.false.,.false.)
      call ga_symmetrize(g_fock)
      end

      subroutine tddft_finalize_fit()
      implicit none
#include "cdft.fh"
#include "errquit.fh"
#include "bas.fh"
#include "global.fh"
      if (cd_bas_han.eq.-99999) return

      if (.not.bas_destroy(cd_bas_han)) then
        call errquit("tddft: cannot destroy cdba",cd_bas_han,BASIS_ERR)
      endif

      if (noio.eq.1) then
        if (.not.ga_destroy(g_2ceri) .or. .not.ga_destroy(g_cdinv))
     $    call errquit("tddft: cannot destroy GAs",0,GA_ERR)   
      endif
      end
