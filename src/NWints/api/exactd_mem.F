c $Id$
*
***********************************************************************
c:tex-% part of the API Internal Routines
c:tex-\subsection{exactd\_mem}
c:tex-This routine computes the memory required by the 
c:tex-McMurchie-Davidson integral derivative code developed at PNNL.
c:tex-This calls specific routines listed below for each integral
c:tex-derivative type.  The data is stored in a common block (apiP.fh)
c:tex-for fast retrieval from the int\_mem routines.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine exactd_mem(rtdb,bases,nbas)
c:tex-\end{verbatim}
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "apiP.fh"
#include "stdio.fh"
#include "bgj.fh"
#include "util.fh"
c
c:: functions
      logical api_is_atomic_basis
      external api_is_atomic_basis
c:: passed
c:tex-\begin{verbatim}
      integer rtdb        ! [input] the RTDB handle
      integer nbas        ! [input] number of basis sets
      integer bases(nbas) ! [input} array of basis set handles
c:tex-\end{verbatim}
c:: local
      integer lmem_2e4c, lmem_1e
      integer lsz_2e4c, lsz_1e
      integer ibas,jbas
      integer nbfmax, nbf
      logical o_do_full
c
c check to see if a full memory check is needed
c
c:tex-
c:tex-{\bf Debugging Note:} using a set directive to set the 
c:tex-variable ``int:full\_mem'' to true will force the more 
c:tex-expensive $O(N^4)$ algorithm to be used to compute the 
c:tex-memory requirements for the 2-electron integral derivatives.
c:tex-
      if (.not.(rtdb_get(rtdb,'int:full_mem',Mt_Log,1,o_do_full)))
     &    o_do_full = .false.
c
      call int_mem_zero()
c
      nbfmax = 0
      do 00100 ibas = 1,nbas
        if (api_is_atomic_basis(bases(ibas))) then
          call int_nbf_max(bases(ibas),nbf)
          nbfmax = max(nbfmax,nbf)
          do 00200 jbas = 1,ibas
            if (api_is_atomic_basis(bases(jbas))) then
c.... nbfmax for 2e3c and 2e4c
              if (bases(jbas).ne.bases(ibas)) then
                call int_nbf_max(bases(ibas),nbf)
                nbfmax = max(nbfmax,nbf)
              endif
              if (o_do_full) then
c.... 2e 4 center full computation
                call emem_d2e4c_full
     &              (bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
              else
c.... 2e 4 center
                call emem_d2e4c
     &              (bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
              endif
              mem_2e4c = max(mem_2e4c, lmem_2e4c)
              isz_2e4c = max(isz_2e4c, lsz_2e4c)
c.... 1e
              call emem_d1e(bases(ibas), bases(jbas), lsz_1e, lmem_1e)
              mem_1e = max(mem_1e, lmem_1e)
              isz_1e = max(isz_1e, lsz_1e)
*debug:          write(luout,*)' ibas|jbas ', ibas, '|', jbas
*debug:          write(luout,*)' mem_2e4c: ',mem_2e4c
*debug:          write(luout,*)' mem_1e  : ',mem_1e
            endif
00200     continue
        endif
00100 continue
c
c.. memory for relativistic one-electron integrals
c
      if (dyall_mod_dir) then
        call emem_d1e_rel(sc_bsh, sc_bsh, lsz_1e, lmem_1e)
        mem_1e = max(mem_1e, lmem_1e)
        isz_1e = max(isz_1e, lsz_1e)
      end if
c
c tmp for ecp derivatives
c
      mem_1e = mem_1e * 3
c
c.. need extra memory for H1 call externalized by int_mem_h1 (Ta)
c
      mem_1e_h1 = mem_1e + nbfmax*nbfmax*3*2
      mem_3ov   = 2*mem_1e_h1  ! WAG
c
c... 2e3c and 2e2c use same scratch space so the maxg is the only thing that changes      
c
      isz_2e3c = 3*3*nbfmax**3
      isz_2e2c = 2*3*nbfmax**2
      mem_2e3c = mem_2e4c
      mem_2e2c = mem_2e4c
c
c for NWChem blocking interface
c
      iszb_2e4c = isz_2e4c  ! one quartet minimum
      memb_2e4c = mem_2e4c + mem_2e4c/5 ! + 20 % 
c
      if (numd_red.gt.numd_tot/33) then
        if (ga_nodeid().eq.0) then
          write(luout,*)
     &        ' integral derivatives have > 33% split blocks'
          write(luout,*)
     &        ' insert the following into your input deck'
          write(luout,*)
     &        ' set "intd:approxmem" integer 9000000'
        endif
      endif
c
      end
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_d1e}
c:tex-This routine determines the maximum buffer and scratch size for the 
c:tex-one electron derivative integrals.
c:tex-This routine should not be called by application code.
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_d1e(brain, ketin, max1e, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "util.fh"
c::functions
      integer int_unint 
      external int_unint
c::external bq functions (mv)
      logical geom_extbq_on
      external geom_extbq_on
      integer  geom_extbq_ncenter
      external  geom_extbq_ncenter

c::passed
c:tex-\begin{verbatim}
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer max1e   ! [output] max size of 1e integrals buffer
      integer memsize ! [output] max size of scratch space for 1e integral evaluation
c:tex-\end{verbatim}
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc, nint
      integer bra_type, ket_type, memdum, bra_geom, ket_geom
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket, nat
      integer bra_fact, ket_fact
      double precision factor
c
#include "bas_exndcf_sfn.fh"
c
      memsize = 0
      max1e = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
      bra_geom  = ibs_geom(bra)
      ket_geom  = ibs_geom(ket)
      if (bra_geom.ne.ket_geom) then
        write(luout,*)'emem_d1e: two different geometries for',
     &         ' derivatives?'
        call errquit('emem_d1e: geom error ',911, GEOM_ERR)
      endif
c
      nat = ncenter(bra_geom)
c
c     add number of external bq's if any (MV)
c     --------------------------------------
      if(geom_extbq_on()) then
        nat = nat + geom_extbq_ncenter() 
      end if
c
c
      do 00100 ibc = 1, bra_ucont
        bra_fact = 1
        ket_fact = 1
c
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        if (bra_gen.gt.1) then
          bra_gen = 1
          bra_fact = 1
        endif
        if (bra_type.lt.0) then
          bra_type = abs(bra_type)
          bra_gen  = 1
        endif
c
        if (bra.eq.ket) then
          ket_lim = ibc
        else
          ket_lim = ket_ucont
        endif
c
        do 00200 ikc = 1, ket_lim
          ket_type = infbs_cont(CONT_TYPE,ikc,ket)
          ket_prim = infbs_cont(CONT_NPRIM,ikc,ket)
          ket_gen  = infbs_cont(CONT_NGEN,ikc,ket)
          if (ket_gen.gt.1) then
            ket_gen = 1
            ket_fact = 1
          endif
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 60 000 000
          nint = int_unint(brain,ibc,ketin,ikc,0,0,0,0)
          call hf1d(
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen,
     &          bra_type, 1, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          ket_prim, ket_gen, ket_type, 2,
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), nat, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          nint, .true., .true., .true., .false., .true.,
     &          dbl_mb(mb_exndcf(1,1)), memdum)
          memsize = max(memsize,(bra_fact*ket_fact*memdum))
00200   continue
00100 continue
c
      call int_nbf_max(brain,max_nbf_bra)
      if (ketin.eq.brain) then
        max_nbf_ket = max_nbf_bra
      else
        call int_nbf_max(ketin,max_nbf_ket)
      endif
c
      max1e = max_nbf_bra*max_nbf_ket*nat*3
      memsize = max(max1e,memsize)
c
      factor = 1.20d00
      memsize = int(dble(memsize)*factor)
      max1e   = int(dble(max1e)*factor)
*debug:      write(luout,'(a,i10)')'emem_d1e: memsize =',memsize
*debug:      write(luout,'(a,i10)')'emem_d1e: max1e   =',max1e
*debug:      write(luout,*)' factor =',factor
c
      end
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_d1e\_rel}
c:tex-This routine determines the maximum buffer and scratch size for the 
c:tex-one electron derivative relativistic integrals.
c:tex-This routine should not be called by application code.
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_d1e_rel(brain, ketin, max1e, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "util.fh"
c::functions
      integer int_unint 
      external int_unint
c::passed
c:tex-\begin{verbatim}
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer max1e   ! [output] max size of 1e integrals buffer
      integer memsize ! [output] max size of scratch space for 1e integral evaluation
c:tex-\end{verbatim}
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc, nint
      integer bra_type, ket_type, memdum, bra_geom, ket_geom
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket, nat
      integer bra_fact, ket_fact
      integer rel_dbg, rel_typ
      double precision factor
c
#include "bas_exndcf_sfn.fh"
c
      memsize = 0
      max1e = 0
c
c     relativistic integral options. rel_typ selects scalar (1),
c     spin-orbit (3) or both (4). Should be passed when spin-orbit 
c     integrals are needed, or included in a common block.
c
      rel_dbg = 0
      rel_typ = 1
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
      bra_geom  = ibs_geom(bra)
      ket_geom  = ibs_geom(ket)
      if (bra_geom.ne.ket_geom) then
        write(luout,*)'emem_d1e_rel: two different geometries for',
     &         ' derivatives?'
        call errquit('emem_d1e_rel: geom error ',911, GEOM_ERR)
      endif
c
      nat = ncenter(bra_geom)
c
      do 00100 ibc = 1, bra_ucont
        bra_fact = 1
        ket_fact = 1
c
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        if (bra_gen.gt.1) then
          bra_gen = 1
          bra_fact = 1
        endif
        if (bra_type.lt.0) then
          bra_type = abs(bra_type)
          bra_gen  = 1
        endif
c
        if (bra.eq.ket) then
          ket_lim = ibc
        else
          ket_lim = ket_ucont
        endif
c
        do 00200 ikc = 1, ket_lim
          ket_type = infbs_cont(CONT_TYPE,ikc,ket)
          ket_prim = infbs_cont(CONT_NPRIM,ikc,ket)
          ket_gen  = infbs_cont(CONT_NGEN,ikc,ket)
          if (ket_gen.gt.1) then
            ket_gen = 1
            ket_fact = 1
          endif
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 20 000 000
          nint = int_unint(brain,ibc,ketin,ikc,0,0,0,0)
          call rel_oneld(
     &        dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &        bra_prim, bra_gen, bra_type, 1, 
     &        dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &        ket_prim, ket_gen, ket_type, 2,
     &        dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)), nat, dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
C .........           doS     doT     doV     canAB
     &        nint, .true., .true., .true., .false.,
C ........     do_nw   do_hnd   nonrel   dryrun
     &        .true., .false., .false., .true.,
     &        dbl_mb(mb_exndcf(1,1)), memdum,rel_dbg,rel_typ)
          memsize = max(memsize,(bra_fact*ket_fact*memdum))
00200   continue
00100 continue
c
      call int_nbf_max(brain,max_nbf_bra)
      if (ketin.eq.brain) then
        max_nbf_ket = max_nbf_bra
      else
        call int_nbf_max(ketin,max_nbf_ket)
      endif
c
      max1e = max_nbf_bra*max_nbf_ket*nat*3
      memsize = max(max1e,memsize)
c
      factor = 1.20d00
      memsize = int(dble(memsize)*factor)
      max1e   = int(dble(max1e)*factor)
*debug:      write(luout,'(a,i10)')'emem_d1e_rel: memsize =',memsize
*debug:      write(luout,'(a,i10)')'emem_d1e_rel: max1e   =',max1e
*debug:      write(luout,*)' factor =',factor
c
      end
c:tex-% part of the API Internal Routines
c:tex-\subsection{emem\_d2e4c}
c:tex-This routine determines the maximum buffer and scratch size for 
c:tex-the given basis sets to compute 2 electron derivative integrals.  
c:tex-The logic used is similar to that of the integral routine.  
c:tex-This routine should not be called by application code.
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_d2e4c(brain, ketin, maxg, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
c::functions
      integer int_unint
      external int_unint
c::passed
c:tex-\begin{verbatim}
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c:tex-\end{verbatim}
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc, nint
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
      integer bra_geom, ket_geom
      integer bra_cent, ket_cent
      integer bra_fact, ket_fact
      logical bra_rel, ket_rel, status_rel
      double precision factor
*debug:-start-debug
*debug:      integer fbra, fket, find
*debug:-start-end
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      memsize = 0
      maxg = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
      bra_geom  = ibs_geom(bra)
      ket_geom  = ibs_geom(ket)
      if (bra_geom.ne.ket_geom) then
        write(luout,*)'emem_d2e4c: two different geometries for',
     &         ' derivatives?'
        call errquit('emem_d2e4c: geom error ',911, GEOM_ERR)
      endif
c 
c relativistic?
c
      status_rel = dyall_mod_dir .and. (.not.nesc_1e_approx)
     &    .and. (brain .eq. ketin) .and. (brain .eq. ao_bsh)
c
      do 00100 ibc = 1, bra_ucont
        bra_fact = 1
        ket_fact = 1
c
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        bra_rel = infbs_cont(CONT_RELLS,ibc,bra) .ne. 0
        bra_cent  = (sf_ibs_cn2ce(ibc,bra))
        if (bra_gen.gt.1) then
          bra_gen = 1
          bra_fact = 1
        endif
        if (bra_type.lt.0) then
          bra_type = abs(bra_type)
          bra_gen  = 1
        endif
c
        if (bra.eq.ket) then
          ket_lim = ibc
        else
          ket_lim = ket_ucont
        endif
c
        do 00200 ikc = 1, ket_lim
          ket_type = infbs_cont(CONT_TYPE,ikc,ket)
          ket_prim = infbs_cont(CONT_NPRIM,ikc,ket)
          ket_gen  = infbs_cont(CONT_NGEN,ikc,ket)
          ket_rel = infbs_cont(CONT_RELLS ,ikc,ket) .ne. 0
          ket_cent  = (sf_ibs_cn2ce(ikc,ket))
          if (ket_gen.gt.1) then
            ket_gen = 1
            ket_fact = 1
          endif
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 20 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulomb
          nint = int_unint(brain,ibc,brain,ibc,ketin,ikc,ketin,ikc)
          if (status_rel .and. (bra_rel .or. ket_rel)) then
            if (ss_one_cent) then
              call rel_2e4cd_sf (
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
c ...........                             canAB   canCD   canPQ
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........     DryRun
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &            bra_rel,ket_rel,ss_one_cent,do_ssss,0)
            else
              call rel_2e4cd_sf (
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,1,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,2,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,3,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,4,
c ...........                             canAB   canCD   canPQ
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........     DryRun
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &            bra_rel,ket_rel,ss_one_cent,do_ssss,0)
            end if
          else
            call hf2d(
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen,
     &          bra_type, 1, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          bra_prim, bra_gen, bra_type, 2,
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), ket_prim, ket_gen,
     &          ket_type, 3, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          ket_prim, ket_gen, ket_type, 4,
     &          dbl_mb(mb_exndcf(1,1)), nint, .false., .false.,
     &          .false.,.true.,dbl_mb(mb_exndcf(1,1)), memdum)
          end if
*debug:-start-debug
*debug:          write(luout,*)
*debug:          write(luout,*)'===============================================',
*debug:     &        '================================='
*debug:          fbra = bra_prim*((bra_type+1)*(bra_type+2)/2)
*debug:          fket = ket_prim*((ket_type+1)*(ket_type+2)/2)
*debug:          find = fbra*fbra*fket*fket
*debug:          write(luout,10000)ibc,ibc,ikc,ikc,
*debug:     &        bra_prim,bra_prim,ket_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:10000     format(1x,'<',4(i4),'>(',4(i3),') nint=',i10,
*debug:     &        ' memory =',i10,' find =',i10,' ratio = ',f10.3)
*debug:-end-debug
          memsize = max(memsize,memdum)
c
          memdum = 20 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . exchange
          nint = int_unint(brain,ibc,brain,ibc,ketin,ikc,ketin,ikc)
          if (status_rel .and. bra_rel .and. ket_rel) then
            if (ss_one_cent) then
              call rel_2e4cd_sf (
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
c ...........                             canAB   canCD   canPQ
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........     DryRun
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &            bra_rel,ket_rel,ss_one_cent,do_ssss,0)
            else
              call rel_2e4cd_sf (
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,1,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,2,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,3,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,4,
c ...........                             canAB   canCD   canPQ
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........     DryRun
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &            bra_rel,ket_rel,ss_one_cent,do_ssss,0)
            end if
          else
            call hf2d(dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen, bra_type,
     &          1, dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), ket_prim, ket_gen, ket_type,
     &          2, dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen, bra_type,
     &          3, dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), ket_prim, ket_gen, ket_type,
     &          4, dbl_mb(mb_exndcf(1,1)), nint,.false.,.false.,
     &          .false., .true., dbl_mb(mb_exndcf(1,1)), memdum)
          end if
*debug:-start-debug
*debug:          find = fbra*fbra*fket*fket
*debug:          write(luout,10000)ibc,ikc,ibc,ikc,
*debug:     &        bra_prim,ket_prim,bra_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:-end-debug
          memsize = max(memsize,memdum)
c
          if (bra.eq.ket) then
            memdum = 20 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bb|bk)
            nint = int_unint(brain,ibc,brain,ibc,brain,ibc,ketin,ikc)
            if (status_rel .and. bra_rel .and. ket_rel) then
              if (ss_one_cent) then
                call rel_2e4cd_sf (
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,bra_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,bra_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,bra_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,ket_cent,
c ...........                             canAB   canCD   canPQ
     &              dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........     DryRun
     &              .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &              bra_rel,ket_rel,ss_one_cent,do_ssss,0)
              else
                call rel_2e4cd_sf (
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,1,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,2,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,3,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,4,
c ...........                             canAB   canCD   canPQ
     &              dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........     DryRun
     &              .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &              bra_rel,ket_rel,ss_one_cent,do_ssss,0)
              end if
            else
              call hf2d(
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen,
     &            bra_type, 1,  dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            bra_prim, bra_gen, bra_type, 2,
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen,
     &            bra_type, 3, dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            ket_prim, ket_gen, ket_type, 4,
     &            dbl_mb(mb_exndcf(1,1)), nint, .false., .false.,
     &            .false., .true., dbl_mb(mb_exndcf(1,1)), memdum)
            end if
*debug:-start-debug
*debug:          find = fbra*fbra*fbra*fket
*debug:          write(luout,10000)ibc,ibc,ibc,ikc,
*debug:     &        bra_prim,bra_prim,bra_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:-end-debug
            memsize = max(memsize,memdum)
c
            memdum = 20 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bk|kk)
            nint = int_unint(brain,ibc,ketin,ikc,ketin,ikc,ketin,ikc)
            if (status_rel .and. bra_rel .and. ket_rel) then
              if (ss_one_cent) then
                call rel_2e4cd_sf (
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,bra_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,ket_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,ket_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,ket_cent,
c ...........                                 canAB   canCD   canPQ
     &              dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........        DryRun
     &              .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &              bra_rel,ket_rel,ss_one_cent,do_ssss,0)
              else
                call rel_2e4cd_sf (
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              bra_prim,bra_gen,bra_type,1,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,2,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,3,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              ket_prim,ket_gen,ket_type,4,
c ...........                                 canAB   canCD   canPQ
     &              dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........        DryRun
     &              .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &              bra_rel,ket_rel,ss_one_cent,do_ssss,0)
              end if
            else
              call hf2d(
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen,
     &            bra_type, 1, dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            ket_prim, ket_gen, ket_type, 2,
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1))
     &            ,dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type
     &            ,3,dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1))
     &            ,dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type
     &            ,4,dbl_mb(mb_exndcf(1,1)),nint,.false.,.false.,
     &            .false.,.true.,dbl_mb(mb_exndcf(1,1)),memdum)
            end if
*debug:-start-debug
*debug:          find = fbra*fket*fket*fket
*debug:          write(luout,10000)ibc,ikc,ikc,ikc,
*debug:     &        bra_prim,ket_prim,ket_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:          write(luout,*)'===============================================',
*debug:     &        '================================='
*debug:-end-debug
            memsize = max(memsize,(bra_fact*ket_fact*memdum))
          endif
00200   continue
00100 continue
c
      call int_nbf_max(brain,max_nbf_bra)
      if (ketin.eq.brain) then
        max_nbf_ket = max_nbf_bra
      else
        call int_nbf_max(ketin,max_nbf_ket)
      endif
      maxg = max_nbf_bra*max_nbf_bra
      maxg = maxg*max_nbf_ket*max_nbf_ket
      maxg = maxg*12
      memsize = max(memsize,maxg)
      factor = 1.50d00
      memsize = int(dble(memsize)*factor)
      maxg    = int(dble(maxg)*factor)
c
*debug:      write(luout,'(a,i10)')'emem_d2e4c: memsize =',memsize
*debug:      write(luout,'(a,i10)')'emem_d2e4c: maxg    =',maxg
*debug:      write(luout,*)' factor =',factor
c
      end
c:tex-% part of the API Internal Routines
c:tex-\subsection{emem\_d2e4c\_full}
c:tex-This routine determines the maximum buffer and scratch size for 
c:tex-the given basis sets to compute 2 electron derivative integrals.  
c:tex-The logic used is similar to that of the integral routine.  
c:tex-This routine should not be called by application code. 
c:tex-this routine computes the complete memory requirements using
c:tex-an $O(N^4$) algorithm
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_d2e4c_full(brain, ketin, maxg, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
c::functions
      integer int_unint
      external int_unint
c::passed
c:tex-\begin{verbatim}
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c:tex-\end{verbatim}
c::local
      integer approx_maxg
      integer approx_memsize
c
      integer bra, bra_ucont, bra_geom, max_nbf_bra
      integer ket, ket_ucont, ket_geom, max_nbf_ket, ket_lim
      integer ibc, i_type, i_prim, i_gen, i_atom
      integer jbc, j_type, j_prim, j_gen, j_atom
      integer kkc, k_type, k_prim, k_gen, k_atom
      integer lkc, l_type, l_prim, l_gen, l_atom
      integer memdum, nint
      logical bra_rel, ket_rel, status_rel
      double precision factor
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      call emem_d2e4c(brain,ketin,approx_maxg,approx_memsize)
      memsize = 0
      maxg = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      write(luout,*)' bra basis ',bs_name(bra)(1:len_bs_name(bra))
      write(luout,*)' ket basis ',bs_name(ket)(1:len_bs_name(ket))
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
      bra_geom  = ibs_geom(bra)
      ket_geom  = ibs_geom(ket)
      if (bra_geom.ne.ket_geom) then
        write(luout,*)'emem_d2e4c: two different geometries for',
     &         ' derivatives?'
        call errquit('emem_d2e4c: geom error ',911, GEOM_ERR)
      endif
c 
c relativistic?
c
      status_rel = dyall_mod_dir .and. (.not.nesc_1e_approx)
     &    .and. (brain .eq. ketin) .and. (brain .eq. ao_bsh)
c
      do ibc = 1, bra_ucont
c
        i_type = infbs_cont(CONT_TYPE ,ibc,bra)
        i_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        i_gen  = infbs_cont(CONT_NGEN ,ibc,bra)
        i_atom = (sf_ibs_cn2ce(ibc,bra))
        if (i_gen.gt.1) then
          i_gen = 1
        endif
        if (i_type.lt.0) then
          i_type = abs(i_type)
          i_gen  = 1
        endif
        do jbc = 1, bra_ucont
          j_type = infbs_cont(CONT_TYPE ,jbc,bra)
          j_prim = infbs_cont(CONT_NPRIM,jbc,bra)
          j_gen  = infbs_cont(CONT_NGEN ,jbc,bra)
          j_atom = (sf_ibs_cn2ce(jbc,bra))
          if (j_gen.gt.1) then
            j_gen = 1
          endif
          if (j_type.lt.0) then
            j_type = abs(j_type)
            j_gen  = 1
          endif
          bra_rel = infbs_cont(CONT_RELLS,ibc,bra) .ne. 0
     &        .and. infbs_cont(CONT_RELLS,jbc,bra) .ne. 0
c
          if (bra.eq.ket) then
            ket_lim = bra_ucont
          else
            ket_lim = ket_ucont
          endif
c
          do kkc = 1, ket_lim
            k_type = infbs_cont(CONT_TYPE ,kkc,ket)
            k_prim = infbs_cont(CONT_NPRIM,kkc,ket)
            k_gen  = infbs_cont(CONT_NGEN ,kkc,ket)
            k_atom = (sf_ibs_cn2ce(kkc,ket))
            if (k_gen.gt.1) then
              k_gen = 1
            endif
            if (k_type.lt.0) then
              k_type = abs(k_type)
              k_gen  = 1
            endif
            do lkc = 1, ket_lim
              l_type = infbs_cont(CONT_TYPE ,lkc,ket)
              l_prim = infbs_cont(CONT_NPRIM,lkc,ket)
              l_gen  = infbs_cont(CONT_NGEN ,lkc,ket)
              l_atom = (sf_ibs_cn2ce(lkc,ket))
              if (l_gen.gt.1) then
                l_gen = 1
              endif
              if (l_type.lt.0) then
                l_type = abs(l_type)
                l_gen  = 1
              endif
              ket_rel = infbs_cont(CONT_RELLS,kkc,ket) .ne. 0
     &            .and. infbs_cont(CONT_RELLS,lkc,ket) .ne. 0
c
              memdum = 20 000 000
              nint = int_unint(brain,ibc,brain,jbc,ketin,kkc,ketin,lkc)
              if (status_rel .and. (bra_rel .or. ket_rel)) then
                call rel_2e4cd_sf (
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              i_prim, i_gen, i_type, i_atom,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              j_prim, j_gen, j_type, j_atom,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              k_prim, k_gen, k_type, k_atom,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              l_prim, l_gen, l_type, l_atom,
c ...........                                 canAB   canCD   canPQ
     &              dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........       DryRun
     &              .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &              bra_rel,ket_rel,ss_one_cent,do_ssss,0)
              else
                call hf2d(
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),
     &            i_prim, i_gen, i_type, i_atom,
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),
     &            j_prim, j_gen, j_type, j_atom,
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),
     &            k_prim, k_gen, k_type, k_atom,
     &            dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),
     &            l_prim, l_gen, l_type, l_atom,
     &            dbl_mb(mb_exndcf(1,1)), nint, .false., .false.,
     &            .false.,.true.,dbl_mb(mb_exndcf(1,1)), memdum)
              end if
              memsize = max(memsize,memdum)
            enddo
          enddo
        enddo
      enddo
c
      call int_nbf_max(brain,max_nbf_bra)
      if (ketin.eq.brain) then
        max_nbf_ket = max_nbf_bra
      else
        call int_nbf_max(ketin,max_nbf_ket)
      endif
      maxg = max_nbf_bra*max_nbf_bra
      maxg = maxg*max_nbf_ket*max_nbf_ket
      maxg = maxg*12
      memsize = max(memsize,maxg)
      factor = dble(maxg)/dble(approx_maxg)
      write(luout,'(a,i10)')  'emem_d2e4c     : maxg    =',approx_maxg
      write(luout,'(a,i10)')  'emem_d2e4c_full: maxg    =',maxg
      write(luout,'(a,f10.3)')'             maxg factor =',factor
      factor = dble(memsize)/dble(approx_memsize)
      write(luout,'(a,i10)')  'emem_d2e4c     : memsize =',
     &    approx_memsize
      write(luout,'(a,i10)')  'emem_d2e4c_full: memsize =',memsize
      write(luout,'(a,f10.3)')'          memsize factor =',factor
      call util_flush(6)
c
      end
