      subroutine dft_getcdt(natoms, nfit, derfit, threshold, compute,
     N     n_batch, n_bmax, iwhat,
     &     n3c_dbl, iwhat_max, i3c_eri,
     T     tvec,g_dens)
      implicit none
#include "errquit.fh"
c
#include "bas.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "schwarz.fh"
#include "msgids.fh"
#include "sym.fh"
      integer natoms
      double precision tvec(*)
      integer nfit
      integer g_dens(nfit)
      double precision threshold
      logical derfit
      logical compute        ! if .true. compute ints for incore
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax         ! max no. batches
      integer iwhat(n_bmax)  ! what batches have been computed [out]
      integer n3c_dbl        ! max no. ints [in]
      integer iwhat_max      !  last batch computed = iwhat(n_batch)
      double precision i3c_ERI(n3c_dbl) ! int values [out]
c
      integer ischw1,nao2
      integer nproc
      double precision fac,Pmat_max,ERI_est,q1
      integer lgvec,igvec,lERI,iERI,lscr,iscr
      integer maxg, mscratch_2e3c,nscr
      integer nintegrals,ninteg_api,ifit
      integer ishc,ishd,ishp
      integer nshbfc,nshbfd,nshbfp
      integer ifirstd,ilastd,ifirstc,ilastc
      integer ifirstp,ilastp
      integer atom_c,atom_d
      integer sh_lo_c,sh_hi_c, sh_lo_d, sh_hi_d, ishd_0, ishd_1
      integer bf_lo_c,bf_hi_c, bf_lo_d, bf_hi_d
      integer nbf_at_c,nbf_at_d
      integer k_p,l_p,icount_d,icount_c
      integer istart,n1
      integer size_atblock
      integer ndone
      integer nxtask,icount
      logical dowork
      double precision dabsmax
      external nxtask,dabsmax
c      
c     Parallelize the loop over the products of AO shells.
c      
      ischw1 = 0
      nproc  = ga_nnodes()
c      

c      
      Nao2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nfit,'gather vector',l_p,k_p))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
      Nao2 = nbf_ao_mxnbf_cn*nbf_ao_mxnbf_cn
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nfit,'gather vector',lgvec,igvec))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
c
      call int_mem_2e3c(maxg, mscratch_2e3c)
      if(.not.MA_Push_Get(MT_Dbl,maxg*nshells_cd,'ERI',lERI,iERI))
     &  call errquit('dft_fitcd: cannot allocate ERI',0, MA_ERR)
      nscr =  mscratch_2e3c 
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0, MA_ERR)
c
c
      n1=0
      if(compute) n_batch=0
      if(.not.incore) n1 = nxtask(nproc,1)
      icount = 0
      ndone = 1
      do atom_c = 1,natoms
        if (.not. bas_ce2cnr( AO_bas_han, atom_c, sh_lo_c, sh_hi_c)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        if (.not. bas_ce2bfr( AO_bas_han, atom_c, bf_lo_c, bf_hi_c)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        nbf_at_c=bf_hi_c - bf_lo_c + 1
c
      do atom_d = 1,atom_c
         if(incore) then
            dowork=iwhat(ndone).eq.((atom_c*(atom_c-1))/2 + atom_d-1)
         else
            dowork=n1.eq.icount
         endif
          if(dowork) then
        if (.not. bas_ce2bfr( AO_bas_han, atom_d, bf_lo_d, bf_hi_d)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        nbf_at_d=bf_hi_d - bf_lo_d + 1
        if(.not.compute) then
c     get DM atom block
        size_atblock= nbf_at_d*nbf_at_c
            do ifit = 1, nfit
               call ga_get(g_dens(ifit), bf_lo_d, bf_hi_d, 
     C              bf_lo_c,bf_hi_c,
     C              dbl_mb(k_p+(ifit-1)*size_atblock), 
     B              (bf_hi_d - bf_lo_d + 1))
            enddo
         endif
        if (.not. bas_ce2cnr( AO_bas_han, atom_d, sh_lo_d, sh_hi_d)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        icount_c=0
        do ishc = sh_lo_c,sh_hi_c
           if (.not. bas_cn2bfr( AO_bas_han, ishc, ifirstc, ilastc)) 
     &          call errquit('Exiting in dft_fitcd',111, BASIS_ERR)
           nshbfc = ilastc-ifirstc+1
           icount_d=0
           do ishd =  sh_lo_d, sh_hi_d
c
            if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, ilastd)) 
     &         call errquit('Exiting in dft_fitcd',1111, BASIS_ERR)
            nshbfd = ilastd-ifirstd+1
            Nao2 = nshbfd*nshbfc
            if (compute) then
               pmat_max=1d0
            else
c
c           Gather the needed density matrix elements.
c
            do ifit = 1, nfit
c            copy pmat into gvec for now
               call dft_copyfrom(
     I              nbf_at_d,nshbfd,nshbfc,
     I              dbl_mb(k_p),dbl_mb(igvec+(ifit-1)*Nao2),
     A              icount_d,icount_c)
            enddo
c            
c           Determine the maximum density matrix element for this product of AO shells.
c
            Pmat_max = dabsmax(Nao2*nfit,DBL_MB(igvec))
         endif 
c
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c
            ERI_est = schwarz_shell(ishc,ishd)
c            
            if( Pmat_max*ERI_est.gt.THRESHOLD )then
c              
              ischw1 = ischw1 + nshells_cd
c             
              istart = iERI

              fac = 2.0d0

              if(atom_c.eq.atom_d) fac = 1.d0
              do ishp = 1,nshells_cd
c                  
c               Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c                  
                if (.not. bas_cn2bfr( CD_bas_han,ishp,ifirstp,ilastp))
     &               call errquit('Exiting in fitcd.',5, BASIS_ERR)
                nshbfp = ilastp-ifirstp+1
c                  
                Nintegrals = nshbfp*Nao2
                Ninteg_api = maxg
c
c               check for use of symmetry
c
c     Must turn symmetry off for derivative fits since derivative
c     densities will not be totally symmetric in general
c
                if (oskel .and. .not.derfit) then
                   if (sym_shell(cd_bas_han, ishp, q1))then
                      call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,
     &                ishd,nscr,Dbl_MB(iscr),Ninteg_api,Dbl_MB(istart))
                      call yscal(nintegrals, q1, Dbl_MB(istart), 1)
                   else
                      call ycopy(nintegrals,0d0,0,Dbl_MB(istart), 1)
                   endif
                else
                   call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,
     &             ishd,nscr,Dbl_MB(iscr),Ninteg_api,Dbl_MB(istart))
                endif                   
c
c
                istart = istart+Nintegrals
c
             enddo
             if(compute) then
c     copy (ahia!) ints to i3c_eri
                if( (ndone+Nao2*nbf_cd-1) .gt. n3c_dbl) then
                   write(LuOut,*) ' nlast n3cdbl ',
     N                  ndone+Nao2*nbf_cd-1, n3c_dbl
                   call errquit('dft_3cincor: blew out of dbls',
     C                  999,UNKNOWN_ERR)
                endif

                call ycopy(Nao2*nbf_cd,Dbl_MB(iERI),1,i3c_eri(ndone))
                ndone = ndone + Nao2*nbf_cd

                n_batch = n_batch + 1
                  if( n_batch.gt.n_bmax)
     &               call errquit('dft_3cincor: out of ints',n_bmax,
     &               UNKNOWN_ERR)
c
                  iwhat(n_batch)=(atom_c*(atom_c-1))/2 + atom_d-1
             else
c               
c             Multiply the matrix of 3-ctr ERIs by the density matrix elements.
c               
              if (nfit.eq.1) then
                 call ygemv('T',Nao2,nbf_cd,fac,Dbl_MB(iERI),Nao2,
     &                Dbl_MB(igvec),1,1.d0,Tvec,1)
              else
                 call errquit('dft_fitcd: not impl',1, BASIS_ERR)
              endif
           endif
            endif
c
c            write(LuOut,*)' Done with: ishc,ishd; n1,n2 ',
c     &                                 ishc,ishd, n1,n2
c            call util_flush(LuOut)
c
c          
        icount_d=icount_d+nshbfd

       enddo ! ishd
        icount_c=icount_c+nshbfc
      enddo ! ishc
      if(.not.incore) n1 = nxtask(nproc,1)
c
          end if
c          
          icount = icount + 1
       enddo ! atom_d
      enddo ! atom_c
  241 continue
c
      if(.not.ma_chop_stack(l_p))
     &  call errquit('dft_fitcd: cannot pop stack',12, MA_ERR)
c
      if(.not.incore) n1 = nxtask(-nproc,1)
c
      if(compute) return
      call ga_mask_sync(incore,.true.)
      call ga_dgop(Msg_gop_Tvec,Tvec,nbf_cd*nfit,'+')
c
c
      if (oskel .and. .not.derfit) then
         do ifit = 1, nfit
            call sym_vec_symmetrize(geom, cd_bas_han,
     &           tvec(1+(ifit-1)*nbf_cd))
         enddo
      endif

      return
      end
