      subroutine cosmo_fock(basis, g)
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "geom.fh"
#include "schwarz.fh"
c     
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     and ADD them into the given global array.
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... requires that ibas = jbas.
c     
c     arguments
c     
      integer basis             !< [Input] bra basis sets
      integer g                !< [Output] GA handle to array 

      integer g_loc
      integer nao2
      integer  maxg_in, mscratch_in
c
      integer  me, nproc, icount,
     &        maxg, mscratch_2e2c, Nintegrals,
     &        ishp, ifirstp, ilastp, nshp,
     &        ishq, ifirstq, ilastq, nshq
      integer lbuf, ibuf, lscr, iscr
      integer geom
      logical odoit
      integer ishq1,ishq2,iat,jat,ishp1,ishp2
      integer iptr,in0,in1,jn0,jn1,nbfa,nbfb
      integer natoms 
      integer vecl
      double precision q2
      double precision threshold
      double precision eri_est
      integer i0,i1,j0,j1,ilo,ihi
      integer iat1,iat2,jat1,jat2
      integer adrc,ldc,j,g_2
      integer mem_lef
      double precision q1,qsym
      integer ishq_2
      integer ijatom
      integer next
      integer nxtask
      external nxtask
      double precision qmax
      logical ocache_save

      call ga_sync()

      ocache_save = ocache
      ocache = .false.

      threshold = 1d-12
      qmax = 1d0

      me = ga_nodeid()
      nproc = ga_nnodes()

      if(.not.ga_duplicate(g,g_loc,'local g')) call
     .     errquit('int_1e_oldga0: dupl failed',0, GA_ERR)
      call ga_zero(g_loc)

      if (.not. bas_geom(basis, geom)) call errquit
     $   ('int_1e_oldga0: basis corrupt?', 0, BASIS_ERR)
c
      call int_mem_1e(maxg_in, mscratch_in)
      mscratch_2e2c = mscratch_in
      maxg = maxg_in

      if (.not. geom_ncent(geom, natoms))
     &     call errquit('ga_get2eri: geom_ncent failed',73, GEOM_ERR)

      if (.not. bas_nbf_ce_max(basis,vecl)) 
     &  call errquit('cosmo_fock: bas_nbf_ce failed',0,BASIS_ERR)
      if (.not.MA_Push_Get(MT_Dbl,mscratch_2e2c,'scr',lscr,iscr))
     &   call errquit('ga_get2eri: cannot allocate scr',2, MA_ERR)
cold      vecl=vecl*vecl*maxg
      vecl=vecl*maxg
      if (.not.MA_Push_Get(MT_Dbl,vecl,'buf',lbuf,ibuf))
     &   call errquit('ga_get2eri: cannot allocate buf',0, MA_ERR)


      ijatom = -1
      next = nxtask(nproc,1)

      do iat=1,natoms
         if (.not. bas_ce2cnr(basis, iat, ishp1, ishp2))
     &        call errquit('cosmo_fock: bas_cn2bfr',0,BASIS_ERR)
         if (.not. bas_ce2bfr(basis, iat, in0, in1))
     &        call errquit('cosmo_fock:  bas_cn2bfr',0,BASIS_ERR)
         nbfa=in1-in0+1
         do jat=1,iat
            if (.not. bas_ce2cnr(basis, jat, ishq1, ishq2))
     &           call errquit('cosmo_fock:bas_cn2bfr',0,BASIS_ERR)
            if (.not. bas_ce2bfr(basis, jat, jn0, jn1))
     &           call errquit('cosmo_fock: bas_cn2bfr',0,BASIS_ERR)
            nbfb=jn1-jn0+1
            ijatom = ijatom + 1
            if (ijatom.eq.next) then

            do ishp=ishp1,ishp2
              if (.not. bas_cn2bfr(basis,ishp,ifirstp,ilastp))
     &            call errquit('ga_get2eri:_cn2bfr',0,BASIS_ERR)
              nshp = ilastp - ifirstp + 1
              call dcopy(nbfb*nshp, 0.0d0,0, dbl_mb(ibuf), 1)
              iptr=ibuf
              do ishq=ishq1,ishq2

                if (.not. bas_cn2bfr(basis,ishq,ifirstq,ilastq))
     &               call errquit('cosmo_fock:bas_cn2bfr',0,BASIS_ERR)
                nshq = ilastq - ifirstq + 1
                nao2 = nshq*nshp

                eri_est = schwarz_shell(ishp,ishq)
                if (eri_est.gt.1.0d-12) then

                Nintegrals=maxg
                call int_1epot(basis,ishq,basis,ishp,mscratch_2e2c,
     &                    Dbl_MB(iscr), Nintegrals,Dbl_MB(iptr))

                if (iat.ne.jat) then
                  call dscal(nao2,2d0,dbl_mb(iptr),1)
                endif

                endif ! schwarz
                iptr=iptr+nao2
              enddo            ! i shell loop
              call ga_put(g_loc,ifirstp,ilastp,jn0,jn1,
     &                    dbl_mb(ibuf),nshp)
c     
            enddo               ! j shell loop
            next = nxtask(nproc,1)
          endif

         enddo                  ! i loop
      enddo                     ! j loop

      call ga_sync()
      next=nxtask(-nproc,1)


      if (.not.ma_chop_stack(lscr))
     $  call errquit('cosmo_fock: could not deallocate scr',0,MA_ERR) 

      call ga_symmetrize(g_loc)

      ocache = ocache_save

      call ga_dadd(1.0d0, g_loc, 1.0d0, g, g)
      call ga_sync()            ! So that no nasty races can result
      if (.not. ga_destroy(g_loc)) call errquit('int_1e_oldga0: gad?',0,
     &     GA_ERR)
      call ga_sync()            ! So that no nasty races can result


      end
