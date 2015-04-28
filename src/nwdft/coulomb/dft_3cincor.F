      subroutine dft_3cincor(n_batch, n_bmax, iwhat, i3c_ERI, 
     &                       n3c_dbl) 
c 
C$Id$
c
      implicit none
#include "errquit.fh"
c    
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches [in]
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl ! max no. ints [in]
      double precision i3c_ERI(n3c_dbl) ! int values [out]
c
#include "global.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "util.fh"
c     
c     Compute the 3-center integrals and store as many of them 
c     as possible in-core.
c
c
c
      double precision THRESHOLD
      integer intdum(1)
      double precision dbldum(1)
c
      THRESHOLD=10.d0**(-itol2e)


      call dcopy(n3c_dbl,0d0,0,i3c_eri, 1)
c      

      call dft_getcdt(ncenters, 1, .false., threshold, .true.,
     N     n_batch, n_bmax, iwhat,
     &     n3c_dbl, i3c_eri,
     T     dbldum, intdum)

      return
      end

      subroutine dft_getcdt(natoms, nfit, derfit, threshold, compute,
     N     n_batch, n_bmax, iwhat,
     &     n3c_dbl, i3c_eri,
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
      double precision i3c_ERI(n3c_dbl) ! int values [out]
c
      integer nao2
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
      integer sh_lo_c,sh_hi_c, sh_lo_d, sh_hi_d
      integer bf_lo_c,bf_hi_c, bf_lo_d, bf_hi_d
      integer nbf_at_c,nbf_at_d
      integer k_p,l_p,icount_d,icount_c
      integer istart,n1
      integer size_atblock
      integer ndone,n_batch_done
      integer nxtask,icount,ipoint
      integer ishd_0,ishd_1
      integer l_at,k_at
      integer atom_c_in,atom_d_in
      logical dowork,use_nxtask,lschw_pmax1,doshell
      double precision dabsmax
      external nxtask,dabsmax
c      
c     Parallelize the loop over the products of AO shells.
c      
      if(derfit) call
     C     errquit('fitcd derfint not impl',5, 0)
      nproc  = ga_nnodes()
c      

c      
      if(.not.compute) then
      Nao2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nfit,'gather vector',l_p,k_p))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
      Nao2 = nbf_ao_mxnbf_cn*nbf_ao_mxnbf_cn
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nfit,'gather vector',lgvec,igvec))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
      endif
c
      call int_mem_2e3c(maxg, mscratch_2e3c)
      if(.not.MA_Push_Get(MT_Dbl,maxg*nshells_cd,'ERI',lERI,iERI))
     &  call errquit('dft_fitcd: cannot allocate ERI',0, MA_ERR)
      nscr =  mscratch_2e3c 
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0, MA_ERR)

      if(.not.MA_Push_Get(MT_int,natoms,'at list',l_at,k_at))
     &  call errquit('dft_3cinc: cannot allocate at list',0, MA_ERR)
c     sort atoms list based on decreasing atomic weight
      call dft_sortatbyz(geom,natoms,int_mb(k_at))
c
c     invalidate iwhat for procs that might remain out of work
      if(compute) then
         do ishp=1,n_bmax
            iwhat(ishp)=-123
         enddo
      endif
c
      n1=0
      use_nxtask=compute.or.(.not.incore)
      n_batch=0
      if(use_nxtask) n1 = nxtask(nproc,1)

      icount = 0
      ndone = 1
      n_batch_done = 1
      ipoint=1
      do atom_c_in = 1,natoms
         atom_c=int_mb(k_at+atom_c_in-1)
        if (.not. bas_ce2cnr( AO_bas_han, atom_c, sh_lo_c, sh_hi_c)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        if (.not. bas_ce2bfr( AO_bas_han, atom_c, bf_lo_c, bf_hi_c)) 
     &       call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        nbf_at_c=bf_hi_c - bf_lo_c + 1
c
      do atom_d_in = 1,atom_c_in
         atom_d=int_mb(k_at+atom_d_in-1)
         if(use_nxtask) then
            dowork=n1.eq.icount
         else
            dowork=iwhat(n_batch_done).eq.
     A           ((atom_c_in*(atom_c_in-1))/2 + atom_d_in-1)
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
     &             call errquit('Exiting in dft_fitcd',1111, BASIS_ERR)
              nshbfd = ilastd-ifirstd+1
c
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c
              ERI_est = schwarz_shell(ishc,ishd)
              if(ERI_est .gt.THRESHOLD) then
              Nao2 = nshbfd*nshbfc
              pmat_max=1d0

              if(.not.compute) then
c
c     Gather the needed density matrix elements.
c
                 do ifit = 1, nfit
c     copy pmat into gvec for now
                    call dft_copyfrom(
     I                   nbf_at_d,nshbfd,nshbfc,
     I                   dbl_mb(k_p),dbl_mb(igvec+(ifit-1)*Nao2),
     A                   icount_d,icount_c)
                 enddo
c     
c     Determine the maximum density matrix element for this product of AO shells.
c     
                 Pmat_max = dabsmax(Nao2*nfit,DBL_MB(igvec))
              endif 
c            
              if( Pmat_max*ERI_est.gt.THRESHOLD )then
c              
                 istart = iERI

                 fac = 2.0d0

                 if(atom_c.eq.atom_d) fac = 1.d0
                 if(.not.incore) then

                    do ishp = 1,nshells_cd
c                  
c     Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c                  
                       if (.not.bas_cn2bfr(CD_bas_han,ishp,
     I                      ifirstp,ilastp))
     &                      call errquit('fitcd.',5, BASIS_ERR)
                       nshbfp = ilastp-ifirstp+1
c                  
                       nintegrals = nshbfp*Nao2
                       Ninteg_api = maxg
c
c     check for use of symmetry
c
c     Must turn symmetry off for derivative fits since derivative
c     densities will not be totally symmetric in general
c
                       q1=1.d0
                       doshell=.true.
                       if(oskel) doshell=sym_shell(cd_bas_han,ishp,q1)

                       if (doshell) then

                          call int_2e3c(CD_bas_han,ishp,AO_bas_han,
     I                         ishc,ishd,nscr,Dbl_MB(iscr),
     N                         Ninteg_api,Dbl_MB(istart))
                          call dscal(nintegrals,q1,dbl_mb(istart),1)
                       else
                          call dcopy(nintegrals,0d0,0,
     D                         dbl_mb(istart),1)
                       endif
                       istart = istart+nintegrals
c     
                    enddo       !ishp

                 endif! end of work for comp. ints
                 if(compute) then
c     incore 3c calculation
c     copy (ahia!) ints to i3c_eri
                    if( (ndone+Nao2*nbf_cd-1) .gt. n3c_dbl) then
                       write(LuOut,*) ' nlast n3cdbl ',
     N                      ndone+Nao2*nbf_cd-1, n3c_dbl
                       call errquit('dft_3cincor: blew out of dbls',
     C                      999,UNKNOWN_ERR)
                    endif

                    call dcopy(Nao2*nbf_cd,dbl_mb(iERI),1,
     I                   i3c_eri(ndone),1)


                 else
c               
c             Multiply the matrix of 3-ctr ERIs by the density matrix elements.
c               
                    if (nfit.eq.1) then
                       if(incore) then

                          call dgemv('T',Nao2,nbf_cd,fac,
     &                         i3c_eri(ipoint),nao2,dbl_mb(igvec),1,
     T                         1.d0,Tvec,1)

                       else
                          call dgemv('T',Nao2,nbf_cd,fac,
     &                         dbl_mb(ieri),nao2,dbl_mb(igvec),1,
     T                         1.d0,Tvec,1)
                       endif
                    else
                       call errquit('dft_fitcd: not impl',1, BASIS_ERR)
                    endif
                 endif
                 ndone  = ndone  + Nao2*nbf_cd
              endif             ! schwartz with real pmax
              ipoint = ipoint + Nao2*nbf_cd
           endif ! schwartz with pmax=1
c          
        icount_d=icount_d+nshbfd

       enddo ! ishd
        icount_c=icount_c+nshbfc
      enddo ! ishc
      if(use_nxtask) then
         n1 = nxtask(nproc,1)
      else
         n_batch_done = n_batch_done + 1
      endif
      if(compute) then
         n_batch = n_batch + 1
         if( n_batch.gt.n_bmax)
     &        call errquit('dft_3cincor: out of ints',n_bmax,
     &        UNKNOWN_ERR)
c     

         iwhat(n_batch)=(atom_c_in*(atom_c_in-1))/2 + atom_d_in-1


      endif
c
          end if
c          
          icount = icount + 1
       enddo ! atom_d
      enddo ! atom_c

c
      if(.not.ma_chop_stack(leri))
     &  call errquit('dft_fitcd: cannot pop stack',12, MA_ERR)

c
      if(use_nxtask) n1 = nxtask(-nproc,1)
      if(compute) return
      if(.not.ma_chop_stack(l_p))
     &  call errquit('dft_fitcd: cannot pop stack',12, MA_ERR)
c
      call ga_mask_sync(.not.use_nxtask,.true.)
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

      subroutine dft_sortatbyz(geom,natoms,klist)
      implicit none
c     sort atoms list based on decreasing atomic weight
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
      integer geom     ! [in]
      integer natoms   ! [in] 
      integer klist(*) ! [out]
      logical geom_cart_get2
      external geom_cart_get2
c
      integer l_xyz,k_xyz,l_z,k_z,l_q,k_q,l_tag,k_tag,found
      integer row,n_rows
      parameter(n_rows=7)
      integer z_value,z_row(n_rows)
#ifdef DEBUB
      integer atom_c,atom_c_in
#endif
      data (z_row(row),row=1,n_rows)/
     1      2,
     2     10,     
     3     18,
     4     36,
     5     54,
     6     86,
     7     103/     

      
      if (.not.MA_Push_Get(MT_int,natoms,'atomic nos',l_z,k_z))
     &   call errquit('xcvdw: cannot allocate atns',0, MA_ERR)
      if (.not.MA_push_get(MT_Dbl,natoms*3,'xyz',l_xyz,k_xyz))
     &   call errquit('xcvdw: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'charge',l_q,k_q))
     &   call errquit('3cinc: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,natoms*16,'tags',l_tag,k_tag))
     &   call errquit('3cinc: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get2(geom, natoms, Byte_MB(k_tag),
     &        Dbl_MB(k_xyz), Dbl_MB(k_q), int_mb(k_z)))
     &   call errquit('3cinc: geom_cart_get failed',0, GEOM_ERR)
      if(.not.ma_chop_stack(l_xyz))
     &  call errquit('dft_fitcd: cannot pop stack',12, MA_ERR)
      found=0
      do row = n_rows,1,-1
         z_value=z_row(row)
         call util_sillysort(z_value,
     R        int_mb(k_z),klist,natoms,found)
      enddo
      call util_sillysort(0,
     R        int_mb(k_z),klist,natoms,found)
      if(natoms.ne.found) then
         write(6,*) ' sorting failed: natoms ',natoms,
     R        ' sorted  ',found
         call errquit('neworder error ',0,0)
      endif
#ifdef DEBUG
      do atom_c_in = 1,natoms
         atom_c=klist(atom_c)
         write(6,'(a,i4,a,i4,a,i4)') ' new list: atom_c_in ',atom_c_in,
     C        ' org indx ',atom_c,' z ',int_mb(k_z+atom_c-1)
      enddo
#endif
      if(.not.ma_chop_stack(l_z))
     &  call errquit('dft_fitcd: cannot pop stack',13, MA_ERR)
      return
      end
