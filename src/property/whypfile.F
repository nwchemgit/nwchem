c*****************************************************************************
      subroutine whypfile(rtdb,
     &                    g_munu_rot,
     &                    g_munu_rot1,
     &                    g_acc2,
     &                    nlst,       
     &                    Ndir_munu,
     &                    Natoms_munu,
     &                    atmnr_munu)
c*****************************************************************************
c
c>>>  Purpose of routine:
c
c>>>  This subroutine reads information from the runtime database
c>>>  and prints it out to a file, name=nbo_fname, which is used
c>>>  as input to the standalone version of the natural bond orbital
c>>>  (NBO) analysis program, gennbo.
c
c>>>  Dependencies/Effects:
c
c     The rtdb is read only and not modified.
c     No common blocks are created.
c     A file, (lfn=25, name=nbo_fname), is created.
c     All variables except for rtdb are local.
c
c>>>  Commentary:
c
c>>>  Eventually, it may be desirable to have the NBO analysis executed
c>>>  as a subroutine of NWchem rather than as a separate stand alone
c>>>  program.  To do this one will have to make a subroutine similar
c>>>  to the one here, but instead of writing the data to an input file
c>>>  it will need to be passed to NBO interface subroutines,  (see the
c>>>  NBO programmers manual for details).  The current strategy was
c>>>  chosen to reduce the binary size of NWchem (the size of the NBO 
c>>>  executable is not negligible), and to take advantage of the fact
c>>>  that gennbo already exists and has been ported to a number of 
c>>>  platforms.
c
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
c#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"

#include "zora.fh" 
c Using g_Ci defined in zora.fh, computed in dft_zora_scale() from dft_zora_utils.F
      integer rtdb  
c
c>>>  handles
c

      integer lbasis, lchg, lcoef,lcoord, leval, lexp, lgeom, lind1, 
     $    lind2, locc, lprim, lptr, lscr, ltags, ltype,
     &    lncomp
      integer g_bo(2), g_fock(2), g_movecs(2),  g_over, g_scr
c
c>>>  indices
      integer ichg, icoef, icoord, ieval, iexp, iind1, iind2, iocc,
     $    iprim, iptr, iscr, itags, itype,
     &    incomp
c
c>>>  basis set objects
      integer atn, charge, maxL, nat, nbf, ncont, ngeno, nmo, 
     $     nprim, sphcart, type
      character*2   sym, elem
      character*255 name, title, trans
c
c>>>  things related to the NBO input file
      integer lfn
      character*255 nbo_fname
c
c>>>  MO vector and calculation-type objects
      character*4   spintype
      character*255 bas_vecs, movecs_in, title_vecs
      character*20 scftype_vecs
      integer  lenocc, nbf_vecs, nset, nmo_vecs(2)
c
c>>>  Functions
      logical int_normalize
      logical  movecs_read, restricted, movecs_read_header
      integer  inp_strlen, ga_create_atom_blocked

      external int_normalize
      external inp_strlen, ga_create_atom_blocked, movecs_read, 
     $     movecs_read_header
c
c>>>  miscellaneous variables
      character*5 cname(5)
      integer i, iat, iat1,icol, icont, 
     &        ilo, ihi, index, ix, iy ,iz, j,
     $        l, len
      logical sshell, pshell, dshell, fshell, gshell
      data cname / 'CS = ', 'CP = ', 'CD = ', 'CF = ', 'CG = ' /
      integer hypopt,ac
      integer l_Ci,k_Ci, ! FA-scaling-zora
     &        l_qscal,k_qscal,noc(2),indq,count_qscal,indx_inc
      integer ipol,noc1,nocc,nelec,nc1
      logical lzora
c     Assuming noc1=noc2 not true if (total spin .ne. 0)
c              (noc_alpha = noc_beta)
      character*11 fname_hyp
      integer nocc_a,nocc_b,idir,jlo,jhi,alo(3),ahi(3),ld(2)
      integer l_buf,k_buf,
     &        l_buf1,k_buf1,
     &        indx   
      integer g_munu_rot,g_munu_rot1,nlst,g_acc2,
     &        Ndir_munu,
     &        Natoms_munu,atmnr_munu(Natoms_munu)
      logical status
      double precision threshold4data,ac_mo
      data threshold4data /1d-40/ ! Threshhold for output data
      if (ga_nodeid().eq.0)
     & write(*,*) 'Using data threshold :',threshold4data
c     --> To store ONLY munu principal components xx,yy,zz
c    Structure of g_munu_rot:
c    n=nbf  Nr. of basis functions
c    val_11_1_1 val_22_1_1 ... val_nn_1_1 val_21_1_1 
c                                         val_31_1_1 val_32_1_1
c                                         val_41_1_1 val_42_1_1 val_43_1_1 
c                                         ...
c    val_11_2_1 val_22_2_1 ... val_nn_2_1 val_21_2_1 
c                                         val_31_2_1 val_32_2_1
c                                         val_41_2_1 val_42_2_1 val_43_2_1 
c                                         ...
c    ...
c    val_11_1_2 val_22_1_2 ... val_nn_1_2 val_21_1_2 
c                                         val_31_1_2 val_32_1_2
c                                         val_41_1_2 val_42_1_2 val_43_1_2
c                                         ...
c    val_11_2_2 val_22_2_2 ... val_nn_2_2 val_21_2_2 
c                                         val_31_2_2 val_32_2_2
c                                         val_41_2_2 val_42_2_2 val_43_2_2
c                                         ...
c    Notation.- val_PQ_R_S
c               PQ, unique munu indices
c               R, direction index = 1,2,3 (in this order) 
c               S, atom number  R=1,nlist
c               Correspondence with EFG eigenvalues:
c               V33 V22 V11
c               -> If multiplied with corresponding 
c                  scaling density matrix it produces
c                  the total EFG eigenvalues in the same order
c                  as NWChem standard output.            
c               nlist total number of atoms selected in
c               keyword of input script: efieldgradZ4 [nlist] ...
c   
                 
c ======= reading from rtdb (nocc_a,nocc_b) === START
c Those are stored in rtdb in dft_inpana.F
         if (.not. rtdb_get(rtdb, 'prop:Nocc_a',mt_int,1,nocc_a))
     $      call errquit('prop_input-EFGZ4-Nocc_a: rtdb_put failed',
     $                   555, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'prop:Nocc_b',mt_int,1,nocc_b))
     $      call errquit('prop_input-EFGZ4-Nocc_b: rtdb_put failed',
     $                   555, RTDB_ERR)
c ======= reading from rtdb (nocc_a,nocc_b) === END
c ------ Lines below are not necessary ------ START
c      lfn = 25
c      if(ga_nodeid().eq.0) then
c        call util_file_name('hyp', .false., .false., nbo_fname)
c        open(lfn,file=nbo_fname,status='unknown')
c        len = inp_strlen(nbo_fname)
c        write(LUout,9250) nbo_fname(1:len)
c      endif 
c ------ Lines below are not necessary ------ END

c ... Fredy: read HYP option key fron RTDB.
      hypopt = 0
      if (.not.rtdb_get(rtdb, 'prop:hypopt', mt_int, 1,hypopt)) 
     &    call errquit('whypfile: hypopt RTDB failed',0, RTDB_ERR)
c ... Fredy end
c
c>>>  load geometry and symmetry info
      if (.not. geom_create(lgeom, 'geometry'))
     $     call errquit('wnbofile: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, lgeom, 'geometry'))
     $     call errquit('wnbofile: no geometry ', 0, RTDB_ERR)
c     
c>>>  load the basis set and get info about it
      if (.not. bas_create(lbasis, 'mo basis'))
     $     call errquit('prop: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, lgeom, lbasis, 'ao basis')) then
         if (.not. bas_rtdb_load(rtdb, lgeom, lbasis, 'mo basis'))
     $        call errquit('wnbofile: no mo or ao basis set', 0,
     &       RTDB_ERR)
      end if
      if (.not. bas_name(lbasis, name, trans))
     $     call errquit('wnbofile: bas_name?', 0, BASIS_ERR)
c
c>>>  Create ga_array for movecs.

      if ( .not. bas_numbf(lbasis,nbf) )
     $     call errquit('wnbofile: bas_numbf failed', 1, BASIS_ERR)
c
c>>>  For now simply set nmo = nbf.

      nmo = nbf
c     
c>>>    Get information on molecular orbital vectors.
c       
      call util_file_name('movecs', .false., .false., movecs_in)
      if(.not.movecs_read_header(movecs_in,title_vecs,bas_vecs,
     $                           scftype_vecs,nbf_vecs,
     &                           nset,nmo_vecs,2))
     $     call errquit('wnbofile:  failed on movecs_read_header',0,
     &       DISK_ERR)

      if(nset.eq.1)then
        restricted=.true.
        spintype = '    '
      else 
        restricted=.false.
        spintype = 'OPEN' 
      endif 
c
c>>>  Extract title, geometry, nuclear charges, and label info.
      if (.not. geom_ncent(lgeom, nat))
     $     call errquit('wnbofile: geom_ncent failed',1, GEOM_ERR)
c       
c>>>  Begin $GENNBO and $NBO keylist data.
c
c ++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++ HYP printout ++++++++++ START
       if (nset.eq.1) then
        nocc=nocc_a
        noc(1)=nocc_a
       else
        nocc = nocc_a+nocc_b
        noc(1)=nocc_a
        noc(2)=nocc_b
       endif
       lzora=.false.
       if (do_zora .and. so_term.eq.0 .and. 
     &     .not.(not_zora_scale)) lzora=.true.
       if (lzora) then  ! spin-free-zora       
        if (.not.ma_alloc_get(mt_dbl,nset*nocc,'Ci',l_Ci,k_Ci)) 
     &     call errquit('whypfile: ma failed',0,MA_ERR)
        if (.not.ma_alloc_get(mt_dbl,nocc,'qscal',l_qscal,k_qscal)) 
     &     call errquit('whypfile: ma failed',0,MA_ERR)
        call ga_get(g_Ci,1,nset,1,nocc,dbl_mb(k_Ci),nset)
        if (nset.eq.1) then
         indx_inc=1
        else if (nset.eq.2) then
         indx_inc=2
        else
         write(*,*) "ERROR in whypfile:: nset=1 or 2 ..."
         stop
        endif
        count_qscal=1
        do i=1,nset
         indx=i ! =1 for alpha (ODD nrs) =2 for beta (even nrs)
         do j=1,noc(i)
           indq=k_qscal+count_qscal-1
           dbl_mb(indq)=dbl_mb(k_Ci+indx-1)
           if (dbl_mb(indq) .ne. 0.0d0)
     &      dbl_mb(indq)=dbl_mb(indq)-1.0d0

c          if (ga_nodeid().eq.0) then
c           write(*,7) i,j,
c     &                dbl_mb(indq),1.0d0/(1.0d0+dbl_mb(indq)),
c     &                1.0d0/(1.0d0+dbl_mb(indq))**0.5d0
c  7        format('whypfile::qscal(',i4,',',i4,')=(',f15.8,',',
c     &             f15.8,',',f15.8,')')
c          endif
          count_qscal=count_qscal+1
          indx=indx+indx_inc
         enddo ! end-loop-j
        enddo ! end-loop-i
       endif ! end-if-lzora
      if ( ga_nodeid().eq.0) then  ! ------- nodeid0----- START
c ... now a loop over the EFG components, idir = 1,3.
c     number of components was given in namelist
c     ==>  nlst=nbf*(nbf+1)
        if (.not.ma_alloc_get(mt_dbl,nlst,'buf',
     &      l_buf,k_buf))
     &      call errquit('whypfile-buf: ma failed',0,MA_ERR)
        nc1=nbf*nocc
        if (.not.ma_alloc_get(mt_dbl,nc1,'buf1',
     &      l_buf1,k_buf1))
     &      call errquit('whypfile-buf: ma failed',0,MA_ERR)
       do iat1=1,Natoms_munu
        iat=atmnr_munu(iat1)
c=====================================
c ==== Generate filename per atom here
c ====================================
        write (fname_hyp,1) iat
 1      format('hypdata.',i0)
        write(*,*) 'Writing filename=',fname_hyp
        open (lfn, file=fname_hyp, status='unknown')
        write (lfn,
     &   '(a,i10,a,i10,a,i10,/a,i10,a,i10,a,i10,a,i10,/a,l3,a/a)') 
     &   '&inputdata nbasis=',nbf, 
     &   ', nocc=',nocc,', nspin=',2,
     &   ', nocc_a=',nocc_a,', nocc_b=',nocc_b,
     &   ', nmo=',nmo, 
     &   ', nelec=',nocc, 
     &   ', threshold = 0.0, lzora=',lzora, 
     &   ', spinorbit= F , label=''HYP'', components=3', 
     &   ', nbofile=''nbodata'', antisym = T, fullmat = F /'
        if (lzora) then  ! spin-free-zora
c     ----move down: g_Ci scaling factors
         write (lfn,'(a)') 'ZORA MO scaling factors'
         write (lfn, '(3e25.16)' ) (dbl_mb(k_qscal+i-1), i = 1,nocc)
        endif
        do idir=1,3
          jlo=nlst*Ndir_munu*(iat1-1)+nlst*(idir-1)+1
          jhi=jlo+nlst-1
          call dcopy(nlst,0.0d0,0,dbl_mb(k_buf),1)
          call ga_get(g_munu_rot,1,1,jlo,jhi,dbl_mb(k_buf),1)
          do i=0,nlst-1
           if (abs(dbl_mb(k_buf+i)) .lt. threshold4data)
     &       dbl_mb(k_buf+i)=0.0d0
          enddo
          write (lfn,'(a,i1)') 'h1mat for direction ',idir
          write (lfn,'(3e25.16)' ) (dbl_mb(k_buf+i), i = 0,nlst-1)
        enddo ! end-loop-directions
c ---- writing perturbed MO + perturbed AO ----- START       
        do idir=1,3
          alo(1)=1
          ahi(1)=nbf
          alo(2)=1
          ahi(2)=nocc
          alo(3)=3*(iat1-1)+idir
          ahi(3)=3*(iat1-1)+idir
          ld(1)=nbf
          ld(2)=nocc
          call dcopy(nc1,0.0d0,0,dbl_mb(k_buf1),1)
          call nga_get(g_acc2,alo,ahi,dbl_mb(k_buf1),ld)
          do i=0,nc1-1
           if (abs(dbl_mb(k_buf1+i)) .lt. threshold4data)
     &       dbl_mb(k_buf1+i)=0.0d0
          enddo
          write (lfn,'(a,i1)') 'c10 for direction ',idir
          write (lfn,'(3e25.16)' ) (dbl_mb(k_buf1+i),i=0,nc1-1)
          jlo=nlst*Ndir_munu*(iat1-1)+nlst*(idir-1)+1
          jhi=jlo+nlst-1
          call dcopy(nlst,0.0d0,0,dbl_mb(k_buf),1)
          call ga_get(g_munu_rot1,1,1,jlo,jhi,dbl_mb(k_buf),1)
          do i=0,nlst-1
           if (abs(dbl_mb(k_buf+i)) .lt. threshold4data)
     &       dbl_mb(k_buf+i)=0.0d0
          enddo
          write (lfn,'(a,i1)') 'h1mat1 for direction ',idir
          write (lfn,'(3e25.16)' ) (dbl_mb(k_buf+i),i=0,nlst-1)
        enddo ! end-loop-directions
c ---- writing perturbed MO + perturbed AO ----- END
        close(lfn) ! close hypdata
       enddo ! end-loop-atoms
         if (.not.ma_free_heap(l_buf)) call
     &      errquit('munu4nbo: ma_free_heap l_buf',0, MA_ERR)
         if (.not.ma_free_heap(l_buf1)) call
     &      errquit('munu4nbo: ma_free_heap l_buf1',0, MA_ERR)
         if (lzora) then
           if (.not.ma_free_heap(l_Ci)) call
     &       errquit('wshldfile: ma_free_heap l_Ci',0, MA_ERR)
           if (.not.ma_free_heap(l_qscal)) call
     &     errquit('wgshiftfile: ma_free_heap l_qscal',0, MA_ERR)
         endif
        endif ! end-if-ga_nodeid()      

        goto 147

       if (ga_nodeid().eq.0) then ! ------ if-ga_nodeid-1--- START
        if(.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
        if(.not.Ma_Push_Get(MT_Dbl, nat*3, 
     &               'coordinates', lcoord, icoord))
     $     call errquit('wnbofile: push_get failed',2, MA_ERR)
        if(.not.Ma_Push_Get(MT_Dbl, nat, 'charges', lchg, ichg))
     $     call errquit('wnbofile: push_get failed',3, MA_ERR)
        if(.not.Ma_Push_Get(MT_Byte, nat*16, 
     &                'center tags', ltags, itags))
     $     call errquit('wnbofile: push_get failed',4, MA_ERR)
        if(.not.geom_cart_get(lgeom, nat, 
     &                 Byte_MB(itags), Dbl_MB(icoord), 
     $     Dbl_MB(ichg)))
     $     call errquit('wnbofile: geom_cart_get failed',
     &                 5, GEOM_ERR)
        if (restricted) then
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'wnbofile: MOs',
     $        0, 0, g_movecs(1))) call errquit('wnbofile: MOs', 0,
     &        GA_ERR)
         call ga_zero(g_movecs(1))
        else
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'wnbofile: alpha MOs',
     $        0, 0, g_movecs(1))) call errquit('wnbofile: alpha MOs', 0,
     &        GA_ERR)
         call ga_zero(g_movecs(1))
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'wnbofile: beta MOs',
     $        0, 0, g_movecs(2))) call errquit('wnbofile: beta MOs', 0,
     &        GA_ERR)
          call ga_zero(g_movecs(2))
        endif
c -------- Creation of NBO auxiliary file -------------- START
        open (lfn, file='nbodata.f47_extra', status='unknown')
c ++++++++++++ EFG printout ++++++++++ END
c ++++++++++++++++++++++++++++++++++++++++
        write(lfn,9000) nat, nbf, spintype
 9000   format('$GENNBO  UPPER  BODM  BOHR  NATOMS=',i5,
     $         '  NBAS=',i5,2x,a4,'  $END')
c
c ... Fredy: here we write more detailed NBO input if requested. 
        if (hypopt.eq.0) then ! default
           write(lfn,'(a)')'$NBO     $END'
        else if (hypopt.eq.1) then ! generate NLMOs, too, and dump info to files
           write(lfn,'(a/a)')'$NBO BNDIDX NBONLMO=W AONBO=W AONLMO=W ',
     $          ' NLMOMO=W STERIC FILE=nbodata $END'
        else if (hypopt.eq.2) then ! generate even more info, incl AO overlap etc.
           write(lfn,'(a/a)')'$NBO BNDIDX NBONLMO=W AONBO=W AONLMO=W ',
     $          ' NLMOMO=W AOMO=W SAO=W STERIC FILE=nbodata $END'
        else ! other options not supported. Fall back to default
           write(lfn,'(a)')'$NBO     $END'
        end if
c ... jochen end
      endif   ! ------ if-ga_nodeid-1--- END
c ------- Destroy ga arrays -----------
c>>>  End $GENNBO and $NBO keylist data.
c
c>>>  Begin $COORD keylist.

      if( ga_nodeid().eq.0)then
        write(lfn,9100)title(1:80)
 9100   format('$COORD',/,a)

        do iat = 1,nat

          ix = 3*(iat-1) + icoord 
          iy = ix+1
          iz = iy+1
          index = itags+(iat-1)*16

          if(.not.geom_tag_to_element(Byte_MB(index),sym,elem,atn))
cedo     $         call errquit('wnbofile: geom_tag_to_element failed',1)
     .         then
          atn=0
          charge=0
       endif

          charge = Dbl_MB(ichg+iat-1)
          write(lfn,9150)atn, charge, Dbl_MB(ix), Dbl_MB(iy), Dbl_MB(iz)

        enddo 

 9150   format(i3,2x,i3,2x,3f20.10)
        write(lfn,9200)
      endif
c
c>>>  End $COORD keylist.
c
c>>>  Begin $BASIS keylist.
c

      if(.not.bas_numcont(lbasis,ncont))
     $     call errquit('wnbofile: bas_numcont failed',0, BASIS_ERR)

      if(.not.MA_Push_Get(MT_INT,ncont,'prim/shell',lprim,iprim))
     $     call errquit('wnbofile: lprim memory alloc. failed',0,
     &       MA_ERR)

      if(.not.MA_Push_Get(MT_INT,ncont+1,'contraction pointer',lptr,
     $     iptr))
     $     call errquit('wnbofile: lptr memory alloc. failed',0,
     &       MA_ERR)

      if(.not.MA_Push_Get(MT_INT,ncont,'contraction type',ltype,
     $     itype))
     $     call errquit('wnbofile: ltype memory alloc. failed',0,
     &       MA_ERR)

      if(.not.MA_Push_Get(MT_INT,ncont,'contraction/shell size',lncomp,
     $     incomp))
     $     call errquit('wnbofile: lncomp memory alloc. failed',0,
     &       MA_ERR)

      if(ga_nodeid().eq.0)then
        write(lfn,'(a)')'$BASIS'
        if(.not.MA_Alloc_Get(MT_Int,nbf,'scratch',lscr,iscr))
     $       call errquit('wnbofile: lscr memory alloc. failed',0,
     &       MA_ERR)
        do i = 1, nmo
          if(.not.bas_bf2ce(lbasis,i,int_mb(iscr+i-1)))
     $         call errquit('wnbofile: bas_bf2ce failed',0, BASIS_ERR)
        enddo 

        write(lfn,'(1x,a)')'CENTER = '
        write(lfn,'(5x,10i6)') (int_mb(iscr+i-1),i=1,nbf)
        
        i = 0
        nprim = 0
        int_mb(iptr) = 1
        maxL = 0
        do icont = 1,ncont
          if(.not.bas_continfo(lbasis,icont,type,int_mb(iprim+icont-1),
     $         ngeno,sphcart))
     $         call errquit('wnbofile: bas_continfo failed',0,
     &       BASIS_ERR)
          INT_mb(itype+icont-1) = type
          maxL = max(maxL,type)
          nprim = nprim + int_mb(iprim+icont-1)
          int_mb(iptr+icont) = 
     $         int_mb(iptr+icont-1) + int_mb(iprim+icont-1)
          index = iscr+i
          if(type.eq.-2)then
c
c>>>        spd shell
c
            int_mb(index)   = 1
            int_mb(index+1) = 101
            int_mb(index+2) = 102
            int_mb(index+3) = 103

            if( sphcart.eq.0)then
              int_mb(index+4) = 201
              int_mb(index+5) = 202
              int_mb(index+6) = 203
              int_mb(index+7) = 204
              int_mb(index+8) = 205
              int_mb(index+9) = 206
              i = i+10
              int_mb(incomp+icont-1) = 10
            elseif(sphcart.eq.1)then
              int_mb(index+4) = 251
              int_mb(index+5) = 252
              int_mb(index+6) = 253
              int_mb(index+7) = 254
              int_mb(index+8) = 255
              i = i+9
              int_mb(incomp+icont-1) = 10
            endif 

          elseif(type.eq.-1) then
c
c>>>        sp shell
c
            int_mb(index)   = 1
            int_mb(index+1) = 101
            int_mb(index+2) = 102
            int_mb(index+3) = 103
            i = i+4
            int_mb(incomp+icont-1) = 4

          elseif(type.eq.0) then
c
c>>>        s shell
c
            int_mb(index)   = 1
            i = i+1
            int_mb(incomp+icont-1) = 1

          elseif(type.eq.1) then
c
c>>>        p shell
c
            do j = 0,2
              int_mb(index+j) = 101+j
            enddo 
            i = i+3
            int_mb(incomp+icont-1) = 3

          elseif(type.eq.2) then
c
c>>>        d shell
c
            if( sphcart.eq.0)then
              do j = 0,5
                int_mb(index+j) = 201+j
              enddo 
              i = i+6
              int_mb(incomp+icont-1) = 6
            elseif(sphcart.eq.1)then

              do j = 0,4
                int_mb(index+j) = 251+j
              enddo 
              i = i+5
              int_mb(incomp+icont-1) = 5
            endif 

          elseif(type.eq.3) then
c
c>>>        f shell
c
            if(sphcart.eq.0)then
              do j = 0,9
                int_mb(index+j) = 301+j
              enddo 
              i = i+10
              int_mb(incomp+icont-1) = 10
            elseif(sphcart.eq.1)then

              do j = 0,6
                int_mb(index+j) = 351+j
              enddo 
              i = i+7
              int_mb(incomp+icont-1) = 7
            endif 

          elseif(type.eq.4) then
c
c>>>        g shell
c
            if(sphcart.eq.0)then
              do j = 0,14
                int_mb(index+j) = 401+j
              enddo 
              i = i+15
              int_mb(incomp+icont-1) = 15
            elseif(sphcart.eq.1)then
              do j = 0,8
                int_mb(index+j) = 451+j
              enddo
              i = i+9
              int_mb(incomp+icont-1) = 9
            endif 

          elseif(type.gt.4) then
            write(6,*)' only up to g functions allowed '
            call errquit('wnbofile: max angular momentum exceeded',0,
     &       BASIS_ERR)
          endif 
          
        enddo 
c
        write(lfn,'(2x,a)')'LABEL = '
        write(lfn,'(5x,10i6)') (int_mb(iscr+i-1),i=1,nbf)
        write(lfn,9200)

        if (.not. ma_free_heap(lscr))
     $       call errquit('wnbofile: lscr free heap failed',0, MA_ERR)
      endif

c
c>>>  End $BASIS keylist.
c
c
c>>> Begin $CONTRACT keylist.

      if (ga_nodeid().eq.0) then
        if(.not.MA_Alloc_Get(MT_dbl,nprim,'exponents',lexp,iexp))
     $       call errquit('wnbofile: lexp memory alloc. failed',0,
     &       MA_ERR)
        if(.not.MA_Alloc_Get(MT_dbl,nprim,'coefs',lcoef,icoef))
     $       call errquit('wnbofile: lcoef memory alloc failed',0,
     &       MA_ERR)

        index = 0
        do icont = 1,ncont
          if(.not.bas_get_exponent(lbasis,icont,dbl_mb(iexp+index)))
     $         call errquit('wnbofile: failed bas_get_exp',0, BASIS_ERR)
          if(.not.bas_get_coeff(lbasis,icont,dbl_mb(icoef+index)))
     $         call errquit('wnbofile: failed bas_get_coeff',0,
     &       BASIS_ERR)
          index = index + int_mb(iprim+icont-1)
        enddo 


        write(lfn,'(a)') '$CONTRACT'
        write(lfn,'(1x,a10,i5)') 'NSHELLS = ', ncont
        write(lfn,'(1x,a10,i5)') '   NEXP = ', nprim
        write(lfn,'(1x,a10)') '  NCOMP = '
        write(lfn,'(5x,10i6)') (int_mb(incomp+i),i=0,(ncont-1))
        write(lfn,'(1x,a10)') '  NPRIM = '
        write(lfn,'(5x,10i6)') (int_mb(iprim+i-1),i=1,ncont)
        write(lfn,'(1x,a10)') '   NPTR = '
        write(lfn,'(5x,10i6)') (int_mb(iptr+i-1),i=1,ncont)
        write(lfn,'(1x,a10)') '    EXP = '
        write(lfn,'(5x,4f20.10)') (dbl_mb(iexp+i-1),i=1,nprim)

        do l = 0,maxL
          write(lfn,'(5x,a)') cname(l+1)
c
c>>>      Reuse array space pointed to by iexp to print out
c>>>      s, p, d, f and g coefficients.
c
          call dcopy(nprim,0.0d0,0,dbl_mb(iexp),1)
          index = 0
          do icont = 1,ncont
            
            type=int_mb(itype+icont-1)

            sshell = (l.eq.0).and.(type.le.0)
            pshell = (l.eq.1).and.((type.eq.1).or.(type.le.-1))
            dshell = (l.eq.2).and.((type.eq.2).or.(type.le.-2))
            fshell = (l.eq.3).and.(type.eq.3)
            gshell = (l.eq.4).and.(type.eq.4)

            ilo = int_mb(iptr+icont-1)
            ihi = int_mb(iptr+icont)-1
c
c>>>        S-type coefficients.
c
            if(sshell)then
              do i = ilo,ihi
                dbl_mb(iexp+i-1)=dbl_mb(icoef+i-1)
              enddo 
            endif 
c
c>>>        P-type coefficients.
c
            if(pshell)then
              do i = ilo,ihi
                dbl_mb(iexp+i-1)=dbl_mb(icoef+i-1)
              enddo 
            endif 
c
c>>>        D-type coefficients.
c
            if(dshell)then
              do i = ilo,ihi
                dbl_mb(iexp+i-1)=dbl_mb(icoef+i-1)
              enddo 
            endif 
c
c>>>        F-type coefficients.
c
            if(fshell)then
              do i = ilo,ihi
                dbl_mb(iexp+i-1)=dbl_mb(icoef+i-1)
              enddo 
            endif 
c
c>>>        G-type coefficients.
c
            if(gshell)then
              do i = ilo,ihi
                dbl_mb(iexp+i-1)=dbl_mb(icoef+i-1)
              enddo 
            endif 

          enddo 
          
          write(lfn,'(5x,4f20.10)') (dbl_mb(iexp+i-1),i=1,nprim)
        enddo 

        write(lfn,9200)

        if (.not. ma_free_heap(lcoef))
     $       call errquit('wnbofile: lcoef free heap failed',0, MA_ERR)

        if (.not. ma_free_heap(lexp))
     $       call errquit('wnbofile: lexp free heap failed',0, MA_ERR)
      endif 
c
c>>> End $CONTRACT keylist.
c
c>>>  Read in molecular orbital vectors and occupation numbers.
c
      lenocc = nset*nmo
      if (.not. ma_push_get(mt_dbl, lenocc, 'wnbofile: mo evals',
     $     leval, ieval)) call errquit
     $     ('wnbofile: insufficient memory?', lenocc, MA_ERR)

      if (.not. ma_push_get(mt_dbl, lenocc, 'wnbofile: mo occ',
     $     locc, iocc)) call errquit
     $     ('wnbofile: mo occ insufficient memory?', lenocc, MA_ERR)

      if (.not. movecs_read(movecs_in, 1, dbl_mb(iocc), dbl_mb(ieval),
     $     g_movecs))
     $     call errquit('scf_movecs_read failed',0, DISK_ERR)

      if(.not.restricted)then
        if (.not. movecs_read(movecs_in, 2,
     $       dbl_mb(iocc+nbf), dbl_mb(ieval+nbf),
     $       g_movecs(2))) then
          call ga_copy(g_movecs(1), g_movecs(2))
          call dcopy(nbf,dbl_mb(iocc),1,dbl_mb(iocc+nbf),1)
          call dcopy(nbf,dbl_mb(ieval),1,dbl_mb(ieval+nbf),1)
        endif
      endif
c
c>>>    Begin $LCAOMO keylist.
c       
      if( ga_nodeid().eq.0)then
        if(.not.MA_Alloc_Get(MT_Dbl,nbf,'scratch',lscr,iscr))
     $       call errquit('wnbofile: lscr memory alloc. failed',0,
     &       MA_ERR)

        write(lfn,'(a)')'$LCAOMO'

        do icol = 1,nbf
          call ga_get(g_movecs(1),1,nmo,icol,icol,Dbl_MB(iscr),1)
          write(lfn,'(5f20.10)')(Dbl_MB(iscr+i),i=0,nbf-1)
        enddo 


        if(.not. restricted) then
          do icol = 1,nbf
            call ga_get(g_movecs(2),1,nmo,icol,icol,Dbl_MB(iscr),1)
            write(lfn,'(5f20.10)')(Dbl_MB(iscr+i),i=0,nbf-1)
          enddo
        endif

        write(lfn,9200)
        if (.not. ma_free_heap(lscr))
     $       call errquit('wnbofile: lscr free heap failed',0, MA_ERR)
      endif
c
c>>>  End LCAOMO keylist.
c
c
c>>>  Create bond order matrix from movecs array for density keylist.
c
c                     nmo   noc       t
c       rho(r,r') = N*SUM   SUM   C  C   chi(r) chi(r')
c                      jk    i     ji ik
c
c                                      /  noc     t  \
c       bond order matrix, gamma   =  | N*SUM C  C    |
c                               ij     \   i   ji ik /
c
      if( restricted) then
*ga:1:0
        if (.not. ga_create(MT_DBL, nmo, nmo, 'wnbofile: bond order', 
     $     0, 0, g_bo(1))) call errquit('wnbofile: g_bo(1)', 0, GA_ERR)
        call ga_zero(g_bo(1))
      else 

*ga:1:0
        if (.not. ga_create(MT_DBL, nmo, nmo, 'wnbofile: a bond order', 
     $     0, 0, g_bo(1))) call errquit('wnbofile: g_bo(1)', 0, GA_ERR)
        call ga_zero(g_bo(1))
        if (.not. ga_create(MT_DBL, nmo, nmo, 'wnbofile: b bond order', 
     $     0, 0, g_bo(2))) call errquit('wnbofile: g_bo(2)', 0, GA_ERR)
        call ga_zero(g_bo(2))

      endif 
c
c>>>  Create the density matrix in the MO representation, 
c>>>  a square matrix with orbital occupation 
c>>>  numbers on the diagonals and 0 elsewhere.  Store the matrix
c>>>  temporarily in the matrix reserved for the bond order matrix.
c

      if(.not.MA_push_get(MT_INT, nmo, 'indices 1', lind1, iind1))
     $     call errquit('wnbofile: lind1',0, MA_ERR)
      if(.not.MA_push_get(MT_INT, nmo, 'indices 2', lind2, iind2))
     $     call errquit('wnbofile: lind2',0, MA_ERR)

      if (ga_nodeid().eq.0)then

        do i=1, nmo
          int_mb( iind1+i-1 ) = i
          int_mb( iind2+i-1 ) = i
        enddo

*        write(6,*)' g_bo B4 zero'
*        call ga_print(g_bo(1))
*        call ga_zero (g_bo(1))
*        write(6,*)' g_bo after zero/B4 scatter'
*        call ga_print(g_bo(1))
        call ga_scatter ( g_bo(1), dbl_mb(iocc), int_mb(iind1),
     $       int_mb(iind2), nmo )

*        write(6,*)' g_bo after scatter'
*        call ga_print(g_bo(1))
      endif

*ga:1:0
      if (.not. ga_create(MT_DBL, nmo, nmo, 'wnbofile: scratch', 
     $     0, 0, g_scr)) call errquit('wnbofile: scratch', 0, GA_ERR)
      call ga_zero(g_scr)
c                                                t
c>>>  Form the product (Movecs)(Gamma_mo)(Movecs)  = Bond order matrix.
c
      call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_movecs(1), 
     $     g_bo(1), 0.0d0, g_scr)
      call ga_dgemm('n', 't', nmo, nmo, nmo, 1.0d0, g_scr, 
     $     g_movecs(1), 0.0d0, g_bo(1))

      if(.not. restricted) then

        if (ga_nodeid().eq.0)then

          do i=1, nmo
            int_mb( iind1+i-1 ) = i
            int_mb( iind2+i-1 ) = i
          enddo

          call ga_scatter ( g_bo(2), dbl_mb(iocc+nmo), int_mb(iind1),
     $         int_mb(iind2), nmo )

        endif

        call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_movecs(2), 
     $       g_bo(2), 0.0d0, g_scr)
        call ga_dgemm('n', 't', nmo, nmo, nmo, 1.0d0, g_scr, 
     $       g_movecs(2), 0.0d0, g_bo(2))
      endif 
c
c>>>  Begin $DENSITY keylist.

      if( ga_nodeid().eq.0)then
        if(.not.MA_Alloc_Get(MT_Dbl,nbf,'scratch',lscr,iscr))
     $       call errquit('wnbofile: lscr memory alloc. failed',0,
     &       MA_ERR)

        write(lfn,'(a)')'$DENSITY'
c
c>>>  Write out lower triangular section of alpha density matrix.
c

        do icol = 1,nbf
          call ga_get(g_bo(1),1,nmo,icol,icol,Dbl_MB(iscr),1)
          write(lfn,'(5f20.10)')(Dbl_MB(iscr+i),i=0,icol-1)
        enddo 


        if(.not. restricted) then
          do icol = 1,nbf
            call ga_get(g_bo(2),1,nmo,icol,icol,Dbl_MB(iscr),1)
            write(lfn,'(5f20.10)')(Dbl_MB(iscr+i),i=0,icol-1)
          enddo
        endif

        write(lfn,9200)
        if (.not. ma_free_heap(lscr))
     $       call errquit('wnbofile: lscr free heap failed',0, MA_ERR)
      endif
c
c>>>  End $DENSITY keylist.
c
c
c>>>  Create overlap matrix.
      call int_init(rtdb,1,lbasis)
      if (.not.int_normalize(rtdb,lbasis)) call errquit
     &    ('wnbofile: int_normalize failed',911, INT_ERR)
      g_over  = ga_create_atom_blocked(lgeom, lbasis,'Overlap')
c
      call ga_zero(g_over)
      call int_1e_ga(lbasis, lbasis, g_over, 'overlap', .false.)
c
c>>>  Begin $OVERLAP keylist.
c
      if( ga_nodeid().eq.0)then
        if(.not.MA_Alloc_Get(MT_Dbl,nbf,'scratch',lscr,iscr))
     $       call errquit('wnbofile: lscr memory alloc. failed',0,
     &       MA_ERR)

        write(lfn,'(a)')'$OVERLAP'
c
c>>>  Write out lower triangular section of overlap matrix.
c
        do icol = 1,nbf
          call ga_get(g_over,1,nmo,icol,icol,Dbl_MB(iscr),1)
          write(lfn,'(5f20.10)')(Dbl_MB(iscr+i),i=0,icol-1)
        enddo 

        write(lfn,9200)
        if (.not. ma_free_heap(lscr))
     $       call errquit('wnbofile: lscr free heap failed',0, MA_ERR)

      endif 
c
c>>>  End $OVERLAP keylist.
c
c
c>>>  Create Fock matrix in AO basis.
c
      if(restricted) then
*ga:1:0
        if (.not. ga_create(MT_DBL, nmo, nmo, 'wnbofile: Fock matrix', 
     $     0, 0, g_fock(1))) call errquit('wnbofile: g_fock(1)', 0,
     &       GA_ERR)
        call ga_zero(g_fock(1))
      else 

*ga:1:0
        if (.not. ga_create(MT_DBL, nmo, nmo, 'wnbofile: alpha Fock', 
     $     0, 0, g_fock(1))) call errquit('wnbofile: g_fock(1)', 0,
     &       GA_ERR)
        call ga_zero(g_fock(1))
*ga:1:0
        if (.not. ga_create(MT_DBL, nmo, nmo, 'wnbofile: beta Fock', 
     $     0, 0, g_fock(2))) call errquit('wnbofile: g_fock(2)', 0,
     &       GA_ERR)
        call ga_zero(g_fock(2))

      endif 
c
c>>>  Create Fock matrix in MO representation, 
c>>>  a square matrix with oribital energies
c>>>  on the diagonals and 0 elsewhere.
c
      if (ga_nodeid().eq.0)then

        do i=1, nmo
          int_mb( iind1+i-1 ) = i
          int_mb( iind2+i-1 ) = i
        enddo

        call ga_scatter ( g_fock(1), dbl_mb(ieval), int_mb(iind1),
     $       int_mb(iind2), nmo )

      endif
c
c>>>  Create Fock Matrix in AO representation.
c                                   t
c>>>  Fock_ao = S*Movecs*Fock_mo*Movecs *S
c
      call ga_dgemm('n', 't', nmo, nmo, nmo, 1.0d0, g_fock(1), 
     $       g_movecs(1), 0.0d0, g_scr)
      call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_movecs(1), 
     $       g_scr, 0.0d0, g_fock(1))
      call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_fock(1), 
     $       g_over, 0.0d0, g_scr)
      call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_over, 
     $       g_scr, 0.0d0, g_fock(1))

      if(.not. restricted) then
        
        if (ga_nodeid().eq.0)then

          do i=1, nmo
            int_mb( iind1+i-1 ) = i
            int_mb( iind2+i-1 ) = i
          enddo

          call ga_scatter ( g_fock(2), dbl_mb(ieval+nmo), int_mb(iind1),
     $         int_mb(iind2), nmo )

        endif

        call ga_dgemm('n', 't', nmo, nmo, nmo, 1.0d0, g_fock(2), 
     $       g_movecs(2), 0.0d0, g_scr)
        call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_movecs(2), 
     $       g_scr, 0.0d0, g_fock(2))
        call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_fock(2), 
     $       g_over, 0.0d0, g_scr)
        call ga_dgemm('n', 'n', nmo, nmo, nmo, 1.0d0, g_over, 
     $       g_scr, 0.0d0, g_fock(2))

      endif 
c
c>>>  Begin $FOCK keylist.

      if( ga_nodeid().eq.0)then
        if(.not.MA_Alloc_Get(MT_Dbl,nbf,'scratch',lscr,iscr))
     $       call errquit('wnbofile: lscr memory alloc. failed',0,
     &       MA_ERR)

        write(lfn,'(a)')'$FOCK'
c
c>>>  Write out lower triangular section of AO Fock matrix.

        do icol = 1,nbf
          call ga_get(g_fock(1),1,nmo,icol,icol,Dbl_MB(iscr),1)
          write(lfn,'(5f20.10)')(Dbl_MB(iscr+i),i=0,icol-1)
        enddo 

        if(.not. restricted) then
          do icol = 1,nbf
            call ga_get(g_fock(2),1,nmo,icol,icol,Dbl_MB(iscr),1)
            write(lfn,'(5f20.10)')(Dbl_MB(iscr+i),i=0,icol-1)
          enddo
        endif

        write(lfn,9200)
        if (.not. ma_free_heap(lscr))
     $       call errquit('wnbofile: lscr free heap failed',0, MA_ERR)
      endif

c
c>>>  End $FOCK keylist.
c
c
c>>>  Print this out for debugging purposes.

      close(lfn)

c -------- Creation of NBO auxiliary file -------------- END
c
c>>>  Free handles and destroy arrays.

      if (restricted) then
        if (.not. ga_destroy(g_bo(1))) 
     $       call errquit('wnbofile: destroy g_bo',0, GA_ERR)
        if (.not. ga_destroy(g_scr)) 
     $       call errquit('wnbofile: destroy g_scr',0, GA_ERR)
        if (.not. ga_destroy(g_fock(1))) 
     $       call errquit('wnbofile: destroy g_fock',0, GA_ERR)
        if (.not. ga_destroy(g_movecs(1))) 
     $       call errquit('wnbofile: destroy g_movecs',0, GA_ERR)
        if (.not. ga_destroy(g_over)) 
     $     call errquit('wnbofile: destroy g_over',0, GA_ERR)
      else
        if (.not. ga_destroy(g_bo(1)))
     $       call errquit('wnbofile: destroy alpha BO matrix', 0,
     &       GA_ERR)
        if (.not. ga_destroy(g_bo(2)))
     $       call errquit('wnbofile: destroy beta BO matrix', 0, GA_ERR)
        if (.not. ga_destroy(g_scr)) 
     $       call errquit('wnbofile: destroy g_scr',0, GA_ERR)
        if (.not. ga_destroy(g_fock(1))) 
     $       call errquit('wnbofile: destroy g_fock',0, GA_ERR)
        if (.not. ga_destroy(g_fock(2))) 
     $       call errquit('wnbofile: destroy g_fock',0, GA_ERR)
        if (.not. ga_destroy(g_movecs(1)))
     $       call errquit('wnbofile: destroy alpha MOs', 0, GA_ERR)
        if (.not. ga_destroy(g_movecs(2)))
     $       call errquit('wnbofile: destroy beta MOs', 0, GA_ERR)
        if (.not. ga_destroy(g_over)) 
     $     call errquit('wnbofile: destroy g_over',0, GA_ERR)
      endif

      if (.not. ma_pop_stack(lind2))
     $     call errquit('wnbofile: lind2 pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(lind1))
     $     call errquit('wnbofile: lind1 pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(locc))
     $     call errquit('wnbofile: locc pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(leval))
     $     call errquit('wnbofile: leval pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(lncomp))
     $     call errquit('wnbofile: lncomp pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(ltype))
     $     call errquit('wnbofile: ltype pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(lptr))
     $     call errquit('wnbofile: lptr pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(lprim))
     $     call errquit('wnbofile: lprim pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(ltags))
     $     call errquit('wnbofile: ltags pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(lchg))
     $     call errquit('wnbofile: lchg pop stack failed',0, MA_ERR)

      if (.not. ma_pop_stack(lcoord))
     $     call errquit('wnbofile: lcoord pop stack failed',0, MA_ERR)
c ------  Deallocate memory ----------

      call int_terminate()

 147  continue

      if (.not. bas_destroy(lbasis)) call errquit
     $     ('wnbofile: basis destroy failed',0, BASIS_ERR)
      if (.not. geom_destroy(lgeom)) call errquit
     $     ('wnbofile: geom destroy failed', 0, GEOM_ERR)
      call ga_sync()

 9200 format('$END')
 9250 format(/,1x,'Input for gennbo program written to file ',a,'.')
      return       
      end
