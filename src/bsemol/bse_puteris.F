      subroutine bse_puteris(erimo,nmo,nstates,nri,nri_me,rilo,ilower,
     $                      iupper,jlower,jupper,eri,blk,polelo,
     $                      polehi,mynpoles)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"

      character(len=2),intent(in) :: blk
      integer,intent(in) :: nmo, nstates, nri, nri_me, rilo
      integer,intent(in) :: ilower, iupper, jlower, jupper
      integer,intent(out) :: eri,polelo,polehi,mynpoles
      double precision,intent(in),dimension(nstates,nmo,nri_me) :: erimo

      character(len=12),parameter :: pname = 'gw_puteris: '
      character(len=6) :: erilabel,idxlabel

      integer ilo,ihi,jlo,jhi,imo,jmo,ipair,npairs,me,rihi
      integer ni,nj,l_tmp,k_tmp,k_loc,ld
      integer l_pair,k_pair,mynpairs,jpair,iproc,nprocs,ulpair,llpair
      integer ibatch,nbatch,myllpair,myulpair,ma_dbl
      integer lmap, kmap

      logical ok
      integer,external :: gw_get_pair

      integer l_oo,k_oo,l_ov,k_ov,l_vv,k_vv
      double precision temp(nri)

      ni = iupper - ilower + 1
      if (ni.lt.1) return

      nj = jupper - jlower + 1
      if (nj.lt.1) return

      rihi = rilo + nri_me - 1
      me = ga_nodeid()
      nprocs = ga_nnodes()

      erilabel = 'eri '//blk
      idxlabel = 'idx '//blk

      npairs = ni*nj
      mynpairs = npairs/nprocs + 1

      if (.not.ga_create(mt_dbl,nri,npairs,erilabel,nri,0,eri))
     $  call errquit(pname//'can''t create '//blk//' block',0,GA_ERR)

      ma_dbl = ma_inquire_stack(mt_dbl)/nri_me
      nbatch = mynpairs/ma_dbl + 1
      mynpairs = npairs/nbatch
      if (mod(npairs,nbatch).ne.0) mynpairs = mynpairs + 1
      mynpairs = min(npairs,mynpairs)

      if (.not.ma_push_get(mt_dbl,mynpairs*nri_me,'temp',l_pair,k_pair))
     $  call errquit(pname//'can''t create temporary array',0,MA_ERR)

      ulpair = 0
      do iproc=0,nprocs-1
        llpair = ulpair + 1
        ulpair = ulpair + npairs/nprocs
        if (iproc.lt.mod(npairs,nprocs)) ulpair = ulpair  + 1
        if (llpair.gt.ulpair) goto 101
        do ibatch=1,nbatch
          myllpair = (ibatch-1)*mynpairs + llpair
          myulpair = min(myllpair + mynpairs - 1,ulpair)
          do ipair=myllpair,myulpair
            jmo = (ipair-1)/ni + 1
            imo = ipair - (jmo-1)*ni + ilower - 1
            jmo = jmo + jlower - 1
            call dcopy(nri_me,erimo(imo,jmo,1),nstates*nmo,
     $                 dbl_mb(k_pair+(ipair-myllpair)*nri_me),1)
          enddo
          call ga_put(eri,rilo,rihi,myllpair,myulpair,dbl_mb(k_pair),
     $                nri_me)
        enddo
 101    continue
        call ga_sync()
      enddo

      ok = ma_chop_stack(l_pair)

      call ga_distribution(eri,me,ld,ld,polelo,polehi)

      if (blk .eq. 'ov') then
        mynpoles = polehi - polelo + 1
      endif

      end subroutine

