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
#include "util.fh"

      integer ndens
      integer nfocks
      integer nder
      integer g_dens
      integer g_fock
      double precision tol2e
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
      integer ityp, ndim, dims(3)

      double precision threshold, az, shmax, cd_coef_max
      double precision pmat_max,eri_est,fac,vcoul_est,q1,qscale,qsym
      double precision maxf, t0

      logical,external :: xc_gotxc
      integer,external :: nxtask_static, nxtask
      double precision, external :: ddot, dabsmax

      ! GA information
      me = ga_nodeid()
      nproc  = ga_nnodes()

      call nga_inquire(g_fock, ityp, ndim, dims)
      nfocks = dims(1)

      ! Whether we read or not CDINV from disk
      iolgc = noio.ne.1

      ! Integral screening threshold
      THRESHOLD = 10.d0**(-itol2e-2)
      THRESHOLD = tol2e

      ! Initialize Fock matrices
      call ga_zero(g_fock)

      ! Skip ERIs and go straight to XC part
      if (sum(abs(jfac(1:ndens))).eq.0d0) goto 100

      ! Memory allocation
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd*ndens,'Tvec',lTvec,iTvec))
     &  call errquit('tddft_fock_fit: cannot allocate Tvec',0, MA_ERR)
      call dcopy(nbf_cd*ndens, 0.D0, 0, Dbl_MB(iTvec), 1)
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd*ndens,'CDCvec',lcdc,icdc))
     &  call errquit('tddft_fock_fit: cannot allocate CDCvec',0, MA_ERR)
      call dcopy(nbf_cd*ndens, 0.D0, 0, Dbl_MB(icdc), 1)

      Nao2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(.not.MA_Push_Get(MT_Dbl,Nao2*ndens,'gather vector',l_p,k_p))
     &  call errquit('tddft_fock_fit: cannot allocate gather vector',0,
     &       MA_ERR)
      Nao2 = nbf_ao_mxnbf_cn*nbf_ao_mxnbf_cn
      if(.not.MA_Push_Get(MT_Dbl,Nao2*ndens,'gathervector',
     &        lgvec,igvec))
     &  call errquit('tddft_fock_fit: cannot allocate gather vector',0,
     &       MA_ERR)

      call int_mem_2e3c(maxg, mscratch_2e3c)
      if(.not.MA_Push_Get(MT_Dbl,maxg*nshells_cd,'ERI',lERI,iERI))
     &  call errquit('dft_fitcd: cannot allocate ERI',0, MA_ERR)
      nscr =  mscratch_2e3c 
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0, MA_ERR)

      ! Limits for NGA_put and NGA_get
      lo(1) = 1
      hi(1) = ndens
      ld(1) = ndens

      ! Part 1. Contract 3c-ERIs with density matrices
      icount = 0
      n1 = nxtask_static(nproc,1)
      do atom_c=1,ncenters
        if (.not. bas_ce2cnr( AO_bas_han, atom_c, sh_lo_c, sh_hi_c)) 
     &       call errquit('Exiting in tddft_fock_fit',110, BASIS_ERR)
        if (.not. bas_ce2bfr( AO_bas_han, atom_c, bf_lo_c, bf_hi_c)) 
     &       call errquit('Exiting in tddft_fock_fit',110, BASIS_ERR)
        nbf_at_c=bf_hi_c - bf_lo_c + 1

        do atom_d=1,atom_c
          dowork=n1.eq.icount

          ! Skip until assigned task
          if(.not.dowork) then
            icount = icount + 1
            cycle
          endif

          if (.not.bas_ce2bfr(AO_bas_han,atom_d,bf_lo_d,bf_hi_d)) 
     &     call errquit('Exiting in tddft_fock_fit',110, BASIS_ERR)
          if (.not.bas_ce2cnr(AO_bas_han,atom_d,sh_lo_d,sh_hi_d)) 
     &     call errquit('Exiting in tddft_fock_fit',110, BASIS_ERR)
          nbf_at_d=bf_hi_d - bf_lo_d + 1

          ! Get DM block for current atom pair
          size_atblock= nbf_at_d*nbf_at_c
          lo(2:3) = (/ bf_lo_d, bf_lo_c /)
          hi(2:3) = (/ bf_hi_d, bf_hi_c /)
          ld(2) = nbf_at_d
          call nga_get(g_dens,lo,hi,dbl_mb(k_p),ld)

          icount_c=0
          do ishc = sh_lo_c,sh_hi_c
             if (.not. bas_cn2bfr(AO_bas_han,ishc,ifirstc,ilastc)) 
     &        call errquit('Exiting in dft_fitcd',111, BASIS_ERR)
             nshbfc = ilastc-ifirstc+1

             icount_d=0
             do ishd =  sh_lo_d, sh_hi_d
c
               if (.not. bas_cn2bfr(AO_bas_han,ishd,ifirstd,ilastd))
     &           call errquit('Exiting in dft_fitcd',1111, BASIS_ERR)
               nshbfd = ilastd-ifirstd+1
c
               ! First schwarz screening
               ERI_est = schwarz_shell(ishc,ishd)
               if(ERI_est.lt.THRESHOLD) then
                 icount_d = icount_d + nshbfd
                 cycle
               endif

               ! Copy DM block for current shell pair
               call tddft_copyfrom(ndens,nbf_at_d,nbf_at_c,nshbfd,
     +          nshbfc,dbl_mb(k_p),dbl_mb(igvec),icount_d,icount_c,
     +          jfac)
c     
               Nao2 = nshbfd*nshbfc

               ! Second Schwarz screening
               Pmat_max = dabsmax(Nao2*ndens,DBL_MB(igvec))
               if( Pmat_max*ERI_est.lt.THRESHOLD )then
                 icount_d = icount_d + nshbfd
                 cycle
               endif

               fac = 2.0d0
               if(atom_c.eq.atom_d) fac = 1.d0

               istart = iERI
               do ishp = 1,nshells_cd
c                  
                 if (.not.bas_cn2bfr(CD_bas_han,ishp,ifirstp,ilastp))
     &             call errquit('fitcd.',5, BASIS_ERR)
                 nshbfp = ilastp-ifirstp+1
c                  
                 nintegrals = nshbfp*Nao2
                 Ninteg_api = maxg
c
                 call int_2e3c(CD_bas_han,ishp,AO_bas_han,
     I                           ishc,ishd,nscr,Dbl_MB(iscr),
     N                           Ninteg_api,Dbl_MB(istart))
                 istart = istart+nintegrals
               enddo

               ! Contract DM blocks with 3c-ERIs
               call dgemm('T','N',nbf_cd,ndens,nao2,fac,
     .               dbl_mb(ieri),nao2,dbl_mb(igvec),nao2,1d0,
     .               dbl_mb(itvec),nbf_cd)

               icount_d = icount_d + nshbfd
            enddo
            icount_c = icount_c  + nshbfc
          enddo
          n1 = nxtask_static(nproc,1)
          icount = icount + 1
        enddo
      enddo
      n1 = nxtask_static(-nproc,1)

c
      ! Collect results from all ranks
      call ga_mask_sync(.false.,.true.)
      call ga_dgop(1,dbl_mb(iTvec),nbf_cd*ndens,'+')
      call ga_sync()

      ! Part 2. Multiply by inverse metric
      ! Read CDINV if necessary
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

      !Multiply T and by CDINV 
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
      endif
      call ga_mask_sync(.true.,.false.)
      call ga_dgop(1,dbl_mb(icdc),nbf_cd*ndens,'+')

      ! Destroy GA if necessary
      if (iolgc) then
        if (.not.ga_destroy(g_cdinv))
     $    call errquit("tddft_fock_fit: destroy GA",g_cdinv,GA_ERR) 
      endif

      ! Obtain maximum values for all Shells
      if(.not.MA_Push_Get(MT_Dbl,nshells_cd,'CD shell max',
     &  lcd_shell,icd_shell))
     &  call errquit('tddft_fock_fit: cannot allocate CD shell max',0,0)
c
      CD_coef_max = 0.D0
      i1 = 1
      do  ishp = 1,nshells_cd
        if( .not. bas_cn2bfr(CD_bas_han,ishp,ishbf,lshbf) )then
          call errquit('tddft_fock_fit',4,0)
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

      ! Part 3. Obtain Fock matrices
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
               maxf = 0d0
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
     I                   ilastd))
     &                   call errquit('Exiting in fitvcoul.',4,0)
                    nshbfd=ilastd-ifirstd+1
                    Nao2 = nshbfc*nshbfd
c           
                    ERI_est=schwarz_shell(ishc,ishd)
                    if( ERI_est*CD_coef_max.lt.THRESHOLD )then
                      icount_d = icount_d + nshbfd
                      cycle
                    endif

                    v_nonzero=.true.
                    qscale=1d0
                    if(ishc.ne.ishd) qscale=2d0
                    istart=ieri
                    nfunz=0
                    call dfill(nao2*ndens,0d0,dbl_mb(igvec),1)
                    do ishp = 1,nshells_cd
c     
                      if( .not.bas_cn2bfr(CD_bas_han, ishp, 
     I                         ifirstp, ilastp))
     &                         call errquit('fitvc',0,BASIS_ERR)
                      nshbf_cd=ilastp-ifirstp+1
                      Vcoul_est=dbl_mb(icd_shell+ishp-1)*ERI_est
                      nintegrals = nao2*nshbf_cd
                      if( Vcoul_est.gt.THRESHOLD )then
                        Ninteg_api=maxg
                        call int_2e3c(CD_bas_han,ishp,ao_bas_han,ishc,
     A                            ishd,nscr,dbl_mb(iscr),ninteg_api,
     &                            Dbl_MB(istart))
                        call dgemm('N','N',nao2,ndens,
     .                    nshbf_cd,qscale,dbl_mb(istart),nao2,
     .                    dbl_mb(icdc+nfunz),nbf_cd,1d0,dbl_mb(igvec),
     .                    nao2)
                      endif ! schw on VCoul_est
c
                      istart = istart + nintegrals
                      nfunz=nfunz+nshbf_cd
c               
                    enddo ! ishp
c
                    call tddft_copyto(ndens,nbf_at_d,nbf_at_c,
     .                    nshbfd,
     .                    nshbfc,dbl_mb(k_p),dbl_mb(igvec),icount_d,
     .                    icount_c)

                    icount_d=icount_d+nshbfd
                  enddo         ! ishd
                  icount_c=icount_c+nshbfc
               enddo            ! ishc
               
               if(v_nonzero) then
                 lo(2:3) = (/ bf_lo_d, bf_lo_c /)
                 hi(2:3) = (/ bf_hi_d, bf_hi_c /)
                 ld(2) = nbf_at_d
                 call nga_put(g_fock,lo,hi,dbl_mb(k_p),ld)
              endif
              n1 = nxtask(nproc,1)
            endif
          icount = icount + 1
        enddo                    ! atom_d
      enddo                     ! atom_c
c
      ! Chop stack
      if(.not.ma_chop_stack(ltvec))
     &  call errquit('dft_fitcd: cannot pop stack',12, MA_ERR)

      n1 = nxtask(-nproc,1)
      call ga_sync()


  100 continue

      if (xc_gotxc().and.nder.eq.0) then
        call xc_newfock(geom,ao_bas_han,jfac,jfac,nbf_ao,
     &        nder, ndens, nfocks, g_fock,g_dens)
      endif

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

      subroutine tddft_copyfrom(nmat,ld1,ld2,row_sub,col_sub,mat,matsub,
     $                          row_ptr,col_ptr,fac)
      implicit none
      integer nmat, ld1, ld2, row_sub, col_sub
      double precision mat(nmat,ld1,ld2), matsub(row_sub,col_sub,nmat)
      integer row_ptr, col_ptr
      double precision fac(*)

      integer imat,j,i

      do imat=1,nmat
        if (fac(imat).eq.0d0) then
          matsub(:,:,imat) = 0d0
          cycle
        endif
      do j=1,col_sub
      do i=1,row_sub
        matsub(i,j,imat) = mat(imat,row_ptr+i,col_ptr+j)
      enddo
      enddo
      enddo

      end

      subroutine tddft_copyto(nmat,ld1,ld2,row_sub,col_sub,mat,matsub,
     $                        row_ptr,col_ptr)
      implicit none
      integer nmat, ld1, ld2, row_sub, col_sub
c      double precision mat(ld1,ld2,nmat), matsub(row_sub,col_sub,nmat)
      double precision mat(nmat,ld1,ld2), matsub(row_sub,col_sub,nmat)
      integer row_ptr, col_ptr

      integer imat,j,i

c      do imat=1,nmat
      do j=1,col_sub
      do i=1,row_sub
      do imat=1,nmat
c        mat(row_ptr+i,col_ptr+j,imat) = matsub(i,j,imat)
        mat(imat,row_ptr+i,col_ptr+j) = matsub(i,j,imat)
      enddo
      enddo
      enddo
c      enddo

      end
