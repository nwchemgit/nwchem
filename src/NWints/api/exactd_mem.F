c $Id: exactd_mem.F,v 1.12 1996-10-20 06:10:04 d3e129 Exp $
*
***********************************************************************
c:tex-% part of the API Internal Routines
c:tex-\subsubsection{exactd\_mem}
c:tex-This routine computes the memory required by the 
c:tex-McMurchie-Davidson integral derivative code developed at PNNL.
c:tex-This calls specific routines listed below for each integral
c:tex-derivative type.  The data is stored in a common block (apiP.fh)
c:tex-for fast retrieval from the int\_mem routines.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine exactd_mem(bases,nbas)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "apiP.fh"
#include "stdio.fh"
c
c:: passed
c:tex-\begin{verbatim}
      integer nbas ! [input] number of basis sets
      integer bases(nbas) ! [input} array of basis set handles
c:tex-\end{verbatim}
c:: local
      integer lmem_2e4c, lmem_1e
      integer lsz_2e4c, lsz_1e
      integer ibas,jbas
      integer nbfmax, nbf
c
      call int_mem_zero()
c
      nbfmax = 0
      do 00100 ibas = 1,nbas
        call int_nbf_max(bases(ibas),nbf)
        nbfmax = max(nbfmax,nbf)
        do 00200 jbas = 1,ibas
c.... nbfmax for 2e3c and 2e4c
          if (bases(jbas).ne.bases(ibas)) then
            call int_nbf_max(bases(ibas),nbf)
            nbfmax = max(nbfmax,nbf)
          endif
c.... 2e 4 center
          call emem_d2e4c(bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
          mem_2e4c = max(mem_2e4c, lmem_2e4c)
          isz_2e4c = max(isz_2e4c, lsz_2e4c)
c.... 1e
          call emem_d1e(bases(ibas), bases(jbas), lsz_1e, lmem_1e)
          mem_1e = max(mem_1e, lmem_1e)
          isz_1e = max(isz_1e, lsz_1e)
*debug:          write(luout,*)' ibas|jbas ', ibas, '|', jbas
*debug:          write(luout,*)' mem_2e4c: ',mem_2e4c
*debug:          write(luout,*)' mem_1e  : ',mem_1e
00200   continue
00100 continue
c
c tmp for ecp derivatives
c
      mem_1e = mem_1e * 3
c
c.. need extra memory for H1 call externalized by int_mem_h1 (Ta)
c
      mem_1e_h1 = mem_1e + nbfmax*nbfmax*3*2
c
c... 2e3c and 2e2c use same scratch space so the maxg is the only thing that changes      
c
      isz_2e3c = nbfmax**3
      isz_2e2c = nbfmax**2
      mem_2e3c = mem_2e4c
      mem_2e2c = mem_2e4c
c
      if (numd_red.gt.numd_tot/33) then
        write(luout,*)' integral derivatives have > 33% split blocks'
        write(luout,*)' insert the following into your input deck'
        write(luout,*)' set "intd:approxmem" integer 9000000'
      endif
c
      end
c:tex-% part of API Internal Routines
c:tex-\subsubsection{emem\_d1e}
c:tex-This routine determines the maximum buffer and scratch size for the 
c:tex-one electron derivative integrals.
c:tex-This routine should not be called by application code.
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_d1e(brain, ketin, max1e, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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
        call errquit('emem_d1e: geom error ',911)
      endif
c
      nat = ncenter(bra_geom)
c
      do 00100 ibc = 1, bra_ucont
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
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
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 15 000 000
          nint = int_unint(brain,ibc,ketin,ikc,0,0,0,0)
          call hf1d(
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen,
     &          bra_type, 1, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          ket_prim, ket_gen, ket_type, 2,
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          nat, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          nint, .true., .true., .true., .false., .true.,
     &          dbl_mb(mb_exndcf(1,1)), memdum)
          memsize = max(memsize,memdum)
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
c
*     write(6,'(a,i10)')'emem_d1e: memsize =',memsize
*     write(6,'(a,i10)')'emem_d1e: max1e   =',max1e
c
      end
c:tex-% part of the API Internal Routines
c:tex-\subsubsection{emem\_d2e4c}
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
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
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
      integer bra_geom, ket_geom, nat
*debug:-start-debug
*debug:      integer fbra, fket, find
*debug:-start-end
c
#include "bas_exndcf_sfn.fh"
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
        call errquit('emem_d2e4c: geom error ',911)
      endif
c
      nat = ncenter(bra_geom)
c
      do 00100 ibc = 1, bra_ucont
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
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
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 15 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
          nint = int_unint(brain,ibc,brain,ibc,ketin,ikc,ketin,ikc)
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
*debug:-start-debug
*debug:          write(luout,*)
*debug:          write(luout,*)'===============================================',
*debug:     &        '================================='
*debug:          fbra = bra_prim*((bra_type+1)*(bra_type+2)/2)
*debug:          fket = ket_prim*((ket_type+1)*(ket_type+2)/2)
*debug:          find = fbra*fbra*fket*fket
*debug:          write(6,10000)ibc,ibc,ikc,ikc,
*debug:     &        bra_prim,bra_prim,ket_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:10000     format(1x,'<',4(i4),'>(',4(i3),') nint=',i10,
*debug:     &        ' memory =',i10,' find =',i10,' ratio = ',f10.3)
*debug:-end-debug
          memsize = max(memsize,memdum)
c
          memdum = 15 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . exchange
          nint = int_unint(brain,ibc,brain,ibc,ketin,ikc,ketin,ikc)
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
*debug:-start-debug
*debug:          find = fbra*fbra*fket*fket
*debug:          write(6,10000)ibc,ikc,ibc,ikc,
*debug:     &        bra_prim,ket_prim,bra_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:-end-debug
          memsize = max(memsize,memdum)
c
          if (bra.eq.ket) then
            memdum = 15 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bb|bk)
            nint = int_unint(brain,ibc,brain,ibc,brain,ibc,ketin,ikc)
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
*debug:-start-debug
*debug:          find = fbra*fbra*fbra*fket
*debug:          write(6,10000)ibc,ibc,ibc,ikc,
*debug:     &        bra_prim,bra_prim,bra_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:-end-debug
            memsize = max(memsize,memdum)
c
            memdum = 15 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bk|kk)
            nint = int_unint(brain,ibc,ketin,ikc,ketin,ikc,ketin,ikc)
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
*debug:-start-debug
*debug:          find = fbra*fket*fket*fket
*debug:          write(6,10000)ibc,ikc,ikc,ikc,
*debug:     &        bra_prim,ket_prim,ket_prim,ket_prim,
*debug:     &        nint,memdum,find,(dble(memdum)/dble(find))
*debug:          write(luout,*)'===============================================',
*debug:     &        '================================='
*debug:-end-debug
            memsize = max(memsize,memdum)
          endif
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_d2e4c:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_d2e4c:bas_nbf_cn_max error ',2)
c
      maxg = max_nbf_bra*max_nbf_bra
      maxg = maxg*max_nbf_ket*max_nbf_ket
      maxg = maxg*12
      
c
*debug:      write(6,'(a,i10)')'emem_d2e4c: memsize =',memsize
*debug:      write(6,'(a,i10)')'emem_d2e4c: maxg    =',maxg
c
      end

