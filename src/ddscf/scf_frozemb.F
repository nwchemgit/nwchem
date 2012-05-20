      subroutine scf_frozemb(rtdb, geom, basis, nbf_ao, g_frozemb)

      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
c
      integer rtdb, geom, basis, nbf_ao, g_frozemb
c
      integer i,j,ncenters
      double precision fock_froz(nbf_ao,nbf_ao),ovlp(nbf_ao,nbf_ao)
      double precision wght_tot

      integer iqsh, istep, nxyz, ncontrset
      integer ixyz, lxyz, icharge, lcharge, itags, ltags
      integer lrqbuf,irqbuf,lqxyz,iqxyz,lqwgt,iqwgt,nqpts,ncube,
     &     ictr_buf, lvemb, ivemb,iqpts

      double precision rad,x,y,z,w,vemb
      integer lbas_cset_info, ibas_cset_info,
     &     lbas_cent_info, ibas_cent_info,
     &     ldocset, idocset,
     &     l_rchi_atom,i_rchi_atom,
     &     l_rq,i_rq,lniz, iniz,
     &     lchi_ao, ichi_ao,
     &     lfock_froz, ifock_froz 
      integer unitno
      character*255 filevemb
      integer inntsize,ddblsize,ok
c
c     preliminaries
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
      ok = 0
      unitno = 77
c
      do i = 1, nbf_ao
        do j = 1, nbf_ao
            ovlp(i,j) = 0.d0  ! overlap helps test the grid
            fock_froz(i,j) = 0.d0 
        enddo
      enddo
c
c     synchronization
      call ga_sync()
c
c     get the number of atoms
      if (.not. geom_ncent(geom, ncenters)) call errquit
     $     ('scf_frozemb: geom_ncent failed', 0, GEOM_ERR)
      if (ncenters .le. 0) return
c
c     read in the embedding potential defined on a grid
      call util_file_name_noprefix('EMBPOT',.false.,.false., filevemb)

      if (ga_nodeid() .eq. 0) then  ! do this on the master node and pack into ga

       write(luout,*) "Frozen embedding file: ",filevemb
c      write(luout,*) "nbf_ao: ",nbf_ao

       open(unitno, file=filevemb)
       read(unitno,*) nqpts
       write(luout,*) "nqpts: ",nqpts
       if (.not.MA_Push_get(mt_dbl,3*nqpts,'qxyz',lqxyz,iqxyz))
     &  call errquit('scf_frozemb: cannot allocate qxyz',0, MA_ERR)
       if (.not.MA_Push_get(mt_dbl,nqpts,'qwgt',lqwgt,iqwgt))
     &  call errquit('scf_frozemb: cannot allocate qwgt',0, MA_ERR)
       if (.not.MA_Push_get(mt_dbl,nqpts,'vemb',lvemb,ivemb))
     &  call errquit('scf_frozemb: cannot allocate vemb',0, MA_ERR)
c
c      map grid points, weights and potential values to local arrays
c       write(6,*) "iqxyz:",iqxyz
c       write(6,*) "iqwgt:",iqwgt
c       write(6,*) "ivemb:",ivemb

       wght_tot = 0.d0
       do  iqpts = 1,nqpts
        read(unitno,*) x,y,z,w,vemb
        dbl_mb(iqxyz+3*(iqpts-1)) = x
        dbl_mb(iqxyz+3*(iqpts-1)+1) = y
        dbl_mb(iqxyz+3*(iqpts-1)+2) = z
        dbl_mb(iqwgt + iqpts-1) = w
        dbl_mb(ivemb + iqpts-1) = vemb
c        !write(luout,*) x,y,z,w,vemb
c        !write(luout,*) dbl_mb(iqxyz+3*(iqpts-1)),
c     !&                 dbl_mb(iqxyz+3*(iqpts-1)+1),
c     !&                 dbl_mb(iqxyz+3*(iqpts-1)+2),
c     !&                 dbl_mb(iqwgt + iqpts-1),
c     !&                 dbl_mb(ivemb + iqpts-1)
        wght_tot = wght_tot + w
       end do
       close(unitno)

c       write(luout,*) "weight total:",wght_tot
c
c      ao basis set info used by xc_eval_basis
       if (.not.bas_numcont(basis, ncontrset))
     &     call errquit('scf_frozemb:bas_numcont',0, BASIS_ERR)
       if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &   lbas_cent_info, ibas_cent_info))
     &   call errquit('scf_frozemb: cannot allocate bas_cent_info',0,
     &       MA_ERR)
       if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &   lbas_cset_info, ibas_cset_info))
     &   call errquit('scf_frozemb: cannot allocate bas_cset_info',0,
     &       MA_ERR)
       call xc_make_basis_info(basis, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
c
       if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('scf_frozemb: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
       do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
       enddo
c
       if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("scf_frozemb:iniz",0, MA_ERR)
       do i= 1, ncenters
         int_mb(iniz+i-1)=1
       enddo
c
       nxyz = 3*ncenters
       if (.not.MA_push_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('scf_frozemb: cannot allocate xyz',0, MA_ERR)
       if (.not.MA_Push_Get(MT_Dbl,ncenters,'charge',lcharge,icharge))
     &   call errquit('scf_frozemb: cannot allocate charge',0, MA_ERR)
       if (.not.MA_Push_Get(MT_Byte,ncenters*16,'tags',ltags,itags))
     &   call errquit('scf_frozemb: cannot allocate tags',0, MA_ERR)
       if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &  call errquit('scf_frozemb: geom_cart_get failed',74, GEOM_ERR)
c
c      compute the basis functions over the grid
          if(.not.MA_Push_get(MT_dbl, ncenters, 'rchi_atom',
     &              l_rchi_atom,i_rchi_atom))
     &              call errquit("scf_frozemb:rchi_atom",0, MA_ERR)
          if(.not.MA_Push_get(MT_dbl, nqpts*ncenters, 'rq',
     &              l_rq,i_rq))
     &              call errquit("scf_frozemb:rq",0, MA_ERR)
          if (.not.MA_Push_Get(mt_dbl, nqpts*nbf_ao,
     &              'chi_ao', lchi_ao, ichi_ao))
     &              call errquit('scf_frozemb: chi_ao',0, MA_ERR)
          if (.not.MA_Push_Get(mt_dbl, nbf_ao*nbf_ao,
     &              'chi_ao', lfock_froz, ifock_froz))
     &              call errquit('scf_frozemb: chi_ao',0, MA_ERR)
c
          call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters)
          call xc_eval_basis(basis, 0, dbl_mb(ichi_ao),
     &              0d0, 0d0, 0d0, dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters,
     &              int_mb(iniz), log_mb(idocset),
     &              int_mb(ibas_cent_info), int_mb(ibas_cset_info))
c          
c         calculate fock and overlap
          call scf_acc_fock(dbl_mb(ichi_ao), dbl_mb(iqwgt), 
     &         nbf_ao, nqpts, dbl_mb(ivemb), dbl_mb(ifock_froz), ovlp)
c
c         pack into a ga and return
          call ga_put(g_frozemb,1,nbf_ao,1,nbf_ao,dbl_mb(ifock_froz),
     &             nbf_ao)
c
c     clear memory
      if(.not.MA_chop_stack(lqxyz))
     &     call errquit("scf_frozemb: chop qxyz", 100, MA_ERR)
c
      end if  ! ga_nodeid
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
      return
      end
c
c     calculate integral dg < mu(g) | V(g) | nu(g) >
      subroutine scf_acc_fock(chi_ao,wght,nbf,npts,pot,fock,ovlp)

      implicit none

#include "stdio.fh"

      integer nbf, npts
      double precision chi_ao(npts,nbf),wght(npts),pot(npts)
      double precision fock(nbf,nbf),ovlp(nbf,nbf)
      integer i, j, k

      do i = 1, nbf
        do j = 1, nbf
          do k = 1, npts
            fock(i,j)=fock(i,j)+chi_ao(k,i)*pot(k)*wght(k)*chi_ao(k,j)
            ovlp(i,j)=ovlp(i,j)+chi_ao(k,i)*wght(k)*chi_ao(k,j)
          enddo
        enddo
      enddo
c
c      do i = 1, nbf
c        do j = 1, nbf
c          if (i.eq.j) write(luout,*) ovlp(i,j),fock(i,j)
c             write(luout,*) i,j,ovlp(i,j)
c        enddo
c      enddo
c      stop
c
      return
      end
c $Id$
