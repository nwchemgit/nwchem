c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Internal routine to compute the memory required by the 
C> McMurchie-Davidson integral code
C>
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{exact\_mem}
c:tex-This routine computes the memory required by the 
c:tex-McMurchie-Davidson integral code developed at PNNL.  
c:tex-This calls specific routines listed below for each 
c:tex-integral type.  The data is stored in a common block 
c:tex-(apiP.fh) for fast retrieval.
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-   
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine exact_mem(rtdb,bases,nbas)
c:tex-\end{verbatim}
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "apiP.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
c
c:: functions
      logical api_is_atomic_basis
      external api_is_atomic_basis
c:: passed
c:tex-\begin{verbatim}
      integer rtdb        !< [Input] the RTDB handle
      integer nbas        !< [Input] number of basis sets
      integer bases(nbas) !< [Input] basis set handles
c:tex-\end{verbatim}
c:: local
      integer lmem_2e4c, lmem_1e
      integer lsz_2e4c, lsz_1e
      integer lmem_2e3c, lsz_2e3c
      integer lmem_2e2c, lsz_2e2c
      integer ibas,jbas, kbas
      integer lsz_3ov,lmem_3ov
      integer mem2emax
      logical o_do_full
c
c
c check to see if a full memory check is needed
c
c:tex-
c:tex-{\bf Debugging Note:} using a set directive to set the 
c:tex-variable ``int:full\_mem'' to true will force the more 
c:tex-expensive $O(N^4)$ algorithm to be used to compute the 
c:tex-memory requirements for the 2-electron integrals.
c:tex-
      if (.not.(rtdb_get(rtdb,'int:full_mem',Mt_Log,1,o_do_full)))
     &    o_do_full = .false.
c
c... ensure that mem_2e{432}c memory reqirements are zero
c
      call int_mem_zero()
c
c     Get maximum memory for 2e integrals from rtdb
c
      if (.not. rtdb_get(rtdb,'int:mem2emax',Mt_Int,1,mem2emax))
     &    mem2emax = 14 000 000

c
      do 00100 ibas = 1,nbas
        if (api_is_atomic_basis(bases(ibas))) then
          do 00200 jbas = 1,ibas
            if (api_is_atomic_basis(bases(jbas))) then
c.... 2e 4 center
              lmem_2e4c = mem2emax
              if (o_do_full) then
                call emem_2e4c_full
     &              (bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
              else
                call emem_2e4c
     &              (bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
              endif
              mem_2e4c = max(mem_2e4c, lmem_2e4c)
              isz_2e4c = max(isz_2e4c, lsz_2e4c)
c.... 2e 3 center
              lmem_2e3c = mem2emax
              call emem_2e3c
     &            (bases(ibas),bases(jbas),lsz_2e3c,lmem_2e3c)
              mem_2e3c = max(mem_2e3c, lmem_2e3c)
              isz_2e3c = max(isz_2e3c, lsz_2e3c)
c.... 2e 2 center
              lmem_2e2c = mem2emax
              call emem_2e2c
     &            (bases(ibas),bases(jbas),lsz_2e2c,lmem_2e2c)
              mem_2e2c = max(mem_2e2c, lmem_2e2c)
              isz_2e2c = max(isz_2e2c, lsz_2e2c)
c.... 1e
              call emem_1e(bases(ibas), bases(jbas), lsz_1e, lmem_1e)
              mem_1e = max(mem_1e, lmem_1e)
              isz_1e = max(isz_1e, lsz_1e)
            endif
00200     continue
        endif
00100 continue
c
c.. memory for relativistic one-electron integrals
c
      if (dyall_mod_dir) then
        call emem_1e_rel(sc_bsh, sc_bsh, lsz_1e, lmem_1e)
        mem_1e = max(mem_1e, lmem_1e)
        isz_1e = max(isz_1e, lsz_1e)
      end if
c
c.. Memory for Douglas-Kroll integrals
c
      if (doug_kroll) then
         call emem_1e_dk(dk_fitbas,dk_fitbas,lsz_1e, lmem_1e)
         mem_1e = max(mem_1e, lmem_1e)
         isz_1e = max(isz_1e, lsz_1e)
      endif
c
c ... jochen:
c     Memory for X2C integrals. exact_2c and x2c_bas are
c     common-block variables in apiP.fh

      if (exact_2c) then
        if (doug_kroll) call errquit(
     &      'exact_mem: x2c and DK not allowed',1,INPUT_ERR)
        call emem_1e_dk(x2c_bas,x2c_bas,lsz_1e, lmem_1e)
         mem_1e = max(mem_1e, lmem_1e)
         isz_1e = max(isz_1e, lsz_1e)
      endif
c
c.. need extra memory for H1 call externalized by int_mem_h1
c
      mem_1e_h1 = mem_1e + isz_1e
c
c... now do loops for 3 center overlap integrals
c
      isz_3ov = 0
      mem_3ov = 0
      do 00300 ibas = 1,nbas
        do 00400 jbas = 1,ibas
          do 00500 kbas = 1,jbas
            call emem_3ov(bases(ibas),bases(jbas),bases(kbas),
     &          lsz_3ov,lmem_3ov)
            isz_3ov = max(isz_3ov,lsz_3ov)
            mem_3ov = max(mem_3ov,lmem_3ov)
00500     continue
00400   continue
00300 continue
c
c new label routines require additional integral block size
c
      mem_2e4c  = mem_2e4c        +       isz_2e4c
      mem_2e3c  = mem_2e3c        +       isz_2e3c
      mem_2e2c  = mem_2e2c        +       isz_2e2c
      mem_1e    = mem_1e          +       isz_1e
      mem_1e_h1 = mem_1e_h1       +       2*isz_1e
      mem_3ov   = mem_3ov         +       isz_3ov
c
c for NWChem blocking interface
c
      iszb_2e4c = isz_2e4c   ! use minimum of one quartet for integrals
      memb_2e4c = mem_2e4c + mem_2e4c/10  ! +10% to be safe
c
#if defined(API_DEBUG)
c.. debug write statements
      write(luout,*) ' mem_2e4c =',mem_2e4c 
      write(luout,*) ' mem_2e3c =',mem_2e3c 
      write(luout,*) ' mem_2e2c =',mem_2e2c 
      write(luout,*) ' mem_1e   =',mem_1e   
      write(luout,*) ' mem_1e_h1=',mem_1e_h1
      write(luout,*) ' mem_3ov  =',mem_3ov  
      write(luout,*)' '
      write(luout,*)' isz_2e4c =',isz_2e4c
      write(luout,*)' isz_2e3c =',isz_2e3c
      write(luout,*)' isz_2e2c =',isz_2e2c
      write(luout,*)' isz_1e   =',isz_1e
      write(luout,*)' isz_3ov  =',isz_3ov
      write(luout,*)' '
#endif
c
      end
C>
C> \brief Internal routine to compute the memory required by the 
C> 3-center overlap integral code
C>
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_3ov}
c:tex-This routine computes the memory for the 3-center overlap 
c:tex-integrals based on the basis sets used.  
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_3ov(ibasin,jbasin,kbasin,lsz_buf,memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "util.fh"
c::function
c::passed
c:tex-\begin{verbatim}
      integer ibasin  !< [Input] basis set handle for ``i'' contractions
      integer jbasin  !< [Input] basis set handle for ``j'' contractions
      integer kbasin  !< [Input] basis set handle for ``k'' contractions
      integer lsz_buf !< [Output] maximum size of integral buffer 
      integer memsize !< [Output] memory needed for scratch array
c:tex-\end{verbatim}
c::local
      integer ibas, jbas, kbas
      integer it_ucont, jt_ucont, kt_ucont
      integer iic, ijc, ikc
      integer i_type, j_type, k_type
      integer i_prim, j_prim, k_prim
      integer memdum, j_lim, k_lim
      integer nbfimax, nbfjmax, nbfkmax
c
#include "bas_exndcf_sfn.fh"
c
      memsize = 0
      nbfimax = 0
      nbfjmax = 0
      nbfkmax = 0
c
      ibas = ibasin + BASIS_HANDLE_OFFSET
      jbas = jbasin + BASIS_HANDLE_OFFSET
      kbas = kbasin + BASIS_HANDLE_OFFSET
c
      it_ucont = infbs_head(HEAD_NCONT,ibas)
      jt_ucont = infbs_head(HEAD_NCONT,jbas)
      kt_ucont = infbs_head(HEAD_NCONT,kbas)
* error checkint
      if (it_ucont.eq.0 .or. jt_ucont.eq.0 .or. kt_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' i handle ',ibasin
        write(luout,*)' i name   ',bs_name(ibas)
        write(luout,*)' i contractions ',it_ucont
        write(luout,*)' j handle ',jbasin
        write(luout,*)' j name   ',bs_name(jbas)
        write(luout,*)' j contractions ',jt_ucont
        write(luout,*)' k handle ',kbasin
        write(luout,*)' k name   ',bs_name(kbas)
        write(luout,*)' k contractions ',kt_ucont
        call errquit('emem_3ov: bad basis info ',911, BASIS_ERR)
      endif
c
      do 00100 iic = 1,it_ucont
        i_type = abs(infbs_cont(CONT_TYPE,iic,ibas))
        i_prim = infbs_cont(CONT_NPRIM,iic,ibas)
c
c
c  ... not SP aware
c
c
c... currently hf3ois does not allow for general contractions
c... so *_gen is not needed
c        i_gen  = infbs_cont(CONT_NGEN,iic,ibas)
c
        j_lim = jt_ucont
        if (ibas.eq.jbas) j_lim = iic
c
        do 00200 ijc = 1, j_lim
          j_type = abs(infbs_cont(CONT_TYPE,ijc,jbas))
          j_prim = infbs_cont(CONT_NPRIM,ijc,jbas)
c**          j_gen  = infbs_cont(CONT_NGEN,ijc,jbas)
c
          k_lim = kt_ucont
          if (kbas.eq.jbas) k_lim = ijc
          do 00300 ikc = 1, k_lim
            k_type = abs(infbs_cont(CONT_TYPE,ikc,kbas))
            k_prim = infbs_cont(CONT_NPRIM,ikc,kbas)
c**            k_gen  = infbs_cont(CONT_NGEN,ikc,kbas)
            
            memdum  = 14 000 000
c
            call hf3OIs(
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),i_prim,i_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),j_prim,j_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),k_prim,k_type,
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum)
c           Jeff: needs a little (150) extra to do higher multipoles with tiny (3-21G) basis sets,
c                 if it still fails just increase 150 by the amount the error message says you are short
            memsize = max(1024,max(memsize,memdum))
00300     continue
00200   continue
00100 continue
c
      call int_nbf_max(ibasin,nbfimax)
      if (ibasin.eq.jbasin) then
        nbfjmax = nbfimax
      else
        call int_nbf_max(jbasin,nbfjmax)
      endif
      if (kbasin.eq.ibasin.or.kbasin.eq.jbasin) then
        if(kbasin.eq.ibasin) then
          nbfkmax = nbfimax
        else
          nbfkmax = nbfjmax
        endif
      else
        call int_nbf_max(kbasin,nbfkmax)
      endif
      lsz_buf = nbfimax*nbfjmax*nbfkmax
      end
C>
C> \brief Internal routine to compute the memory required by the 
C> 1-electron integral codes
C>
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_1e}
c:tex-This routine computes the memory for any type of 1e 
c:tex-integrals based on the basis sets used.  
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_1e(brain, ketin, max1e, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "util.fh"
c::passed
c:tex-\begin{verbatim}
      integer brain   !< [Input] bra basis set handle
      integer ketin   !< [Input] ket basis set handle
      integer max1e   !< [Output] max size of 1e integrals buffer
      integer memsize !< [Output] max size of scratch space for 1e integral evaluation
c:tex-\end{verbatim}
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
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
*error checking
      if (bra_ucont.eq.0 .or. ket_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' bra handle ',brain
        write(luout,*)' bra name   ',bs_name(bra)
        write(luout,*)' bra contractions ',bra_ucont
        write(luout,*)' ket handle ',ketin
        write(luout,*)' ket name   ',bs_name(ket)
        write(luout,*)' ket contractions ',ket_ucont
        call errquit('emem_1e: bad basis info ',911, BASIS_ERR)
      endif
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
          memdum = 14 000 000
*
          call hf1(dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),1,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),1,.true.,.true.,.true.,.false.,
     &          .true.,dbl_mb(mb_exndcf(1,1)),memdum)
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
      max1e = max_nbf_bra*max_nbf_ket
c... insure scratch is large enough to hold 3*max1e
      memsize = max(3*max1e,memsize)
c... make it at least 1500
      memsize = max(1500,memsize)
c
#if defined(API_DEBUG)
      write(luout,'(a,i10)')'emem_1e: memsize =',memsize
      write(luout,'(a,i10)')'emem_1e: max1e   =',max1e
#endif
c
      end
C>
C> \brief Internal routine to compute the memory required by the 
C> Douglas-Kroll integral codes for the fitting basis
C>
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_1e\_dk}
c:tex-This routine computes the memory for the requested type of  
c:tex-Douglas-Kroll integrals for the given fitting basis sets.
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_1e_dk(brain, ketin, max1e, memsize)
c:tex-\end{verbatim}
      implicit none
#include "apiP.fh"
#include "util.fh"
c::passed
c:tex-\begin{verbatim}
      integer brain   !< [Input] bra basis set handle
      integer ketin   !< [Input] ket basis set handle
      integer max1e   !< [Output] max size of 1e integrals buffer
      integer memsize !< [Output] max size of scratch space for 1e integral evaluation
c:tex-\end{verbatim}
c::local
      integer lsz_1e, lmem_1e
c
      max1e = 0
      memsize = 0
      call emem_1e(brain, ketin, lsz_1e, lmem_1e)
      max1e = max(max1e, lmem_1e)
      memsize = max(memsize, lsz_1e)
      call emem_1e_pvp(brain, ketin, lsz_1e, lmem_1e, 1)
      max1e = max(max1e, lmem_1e)
      memsize = max(memsize, lsz_1e)
      if (dk_form_type .eq. 3) then
         call emem_1e_pvp(brain, ketin, lsz_1e, lmem_1e, 3)
         max1e = max(max1e, lmem_1e)
         memsize = max(memsize, lsz_1e)
      endif
c
      end
c
C>
C> \brief Internal routine to compute the memory required by the 
C> relativistic integral codes
C>
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_1e\_rel}
c:tex-This routine computes the memory for any type of relativistic 1e 
c:tex-integrals based on the basis sets used.  
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_1e_rel(brain, ketin, max1e, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "rel_nwc.fh"
#include "util.fh"
c::passed
c:tex-\begin{verbatim}
      integer brain   !< [Input] bra basis set handle
      integer ketin   !< [Input] ket basis set handle
      integer max1e   !< [Output] max size of 1e integrals buffer
      integer memsize !< [Output] max size of scratch space for 1e integral evaluation
c:tex-\end{verbatim}
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
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
*error checking
      if (bra_ucont.eq.0 .or. ket_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' bra handle ',brain
        write(luout,*)' bra name   ',bs_name(bra)
        write(luout,*)' bra contractions ',bra_ucont
        write(luout,*)' ket handle ',ketin
        write(luout,*)' ket name   ',bs_name(ket)
        write(luout,*)' ket contractions ',ket_ucont
        call errquit('emem_1e_rel: bad basis info ',911, BASIS_ERR)
      endif
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
        do 00200 ikc = 1, ket_ucont
          ket_type = infbs_cont(CONT_TYPE,ikc,ket)
          ket_prim = infbs_cont(CONT_NPRIM,ikc,ket)
          ket_gen  = infbs_cont(CONT_NGEN,ikc,ket)
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 14 000 000
*
          call rel_onel(dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &        bra_prim,bra_gen,bra_type,
     &        dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &        ket_prim,ket_gen,ket_type,
     &        dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)),1,
     &        dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &        dbl_mb(mb_exndcf(1,1)),1,.true.,.true.,.true.,.false.,
     &        .true.,.false.,.false.,.true.,dbl_mb(mb_exndcf(1,1)),
     &        memdum,rel_dbg,rel_typ)
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
      max1e = max_nbf_bra*max_nbf_ket
c... insure scratch is large enough to hold 3*max1e
      memsize = max(3*max1e,memsize)
c
#if defined(API_DEBUG)
      write(luout,'(a,i10)')'emem_1e_rel: memsize =',memsize
      write(luout,'(a,i10)')'emem_1e_rel: max1e   =',max1e
#endif
c
      end
C>
C> \brief Internal routine to compute the memory required by the 
C> 2-electron 2-center integral code
C>
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_2e4c}
c:tex-This routine computes the memory required by the McMurchie-Davidson 
c:tex-algorithm for 4 center two electron integrals based on the 
c:tex-basis sets used.  
c:tex-
c:tex-The exact algorithm is an N$^4$ and an N$^2$ approximate algorithm 
c:tex-is used here.  Exact memory is computed for integrals over the 
c:tex-following classes:
c:tex-\begin{itemize}
c:tex-\item Coulumb $(ii|jj)$
c:tex-\item Exchange $(ij|ij)$
c:tex-\item Triplet $(ii|ij)$
c:tex-\item Triplet $(ij|jj)$
c:tex-\end{itemize}
c:tex-An additional 10\% is added to the maximum exact memory 
c:tex-computation for each of these classes.  Additional classes that 
c:tex-have been periodically problematic are:
c:tex-\begin{itemize}
c:tex-\item $(ij|jk)$
c:tex-\item $(ij|kk)$
c:tex-\end{itemize}
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_2e4c(brain, ketin, maxg, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
c::functions
      logical cando_sp_u
      external cando_sp_u
c::passed
c:tex-\begin{verbatim}
      integer brain   !< [Input] bra basis set handle
      integer ketin   !< [Input] ket basis set handle
      integer maxg    !< [Output] max size of 2e integrals buffer
      integer memsize !< [Output] max size of scratch space for 2e integral evaluation
c:tex-\end{verbatim}
c     NB. On input, memsize contains the (dummy) memory limit for 2e integrals
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum, memmax
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen, bra_cent, ket_cent
      integer max_nbf_bra, max_nbf_ket
c
      logical bra_cando_sp, bra_rel
      logical ket_cando_sp, ket_rel
      logical status_rel
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      memmax = memsize
      memsize = 0
      maxg = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
*error checking
      if (bra_ucont.eq.0 .or. ket_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' bra handle ',brain
        write(luout,*)' bra name   ',bs_name(bra)
        write(luout,*)' bra contractions ',bra_ucont
        write(luout,*)' ket handle ',ketin
        write(luout,*)' ket name   ',bs_name(ket)
        write(luout,*)' ket contractions ',ket_ucont
        call errquit('emem_2e4c: bad basis info ',911, BASIS_ERR)
      endif
c 
c relativistic?
c
      status_rel = dyall_mod_dir .and. (.not.nesc_1e_approx)
     &    .and. (brain .eq. ketin) .and. (brain .eq. ao_bsh)
c
      do 00100 ibc = 1, bra_ucont
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        bra_rel = infbs_cont(CONT_RELLS ,ibc,bra) .ne. 0
        bra_cent  = (sf_ibs_cn2ce(ibc,bra))
c
        if (bra_type.lt.0) then
          bra_type = abs(bra_type)
          bra_gen  = 1
        endif
        bra_cando_sp = cando_sp_u(brain,ibc,0)
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
c
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
          ket_cando_sp = cando_sp_u(ketin,ikc,0)
          if (bra_cando_sp.and.ket_cando_sp) goto 00201
c
          memdum = memmax
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulomb
          if (status_rel .and. (bra_rel .or. ket_rel)) then
            call rel_2e4c_sf (
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          bra_prim,bra_gen,bra_type,bra_cent,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          bra_prim,bra_gen,bra_type,bra_cent,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          ket_prim,ket_gen,ket_type,ket_cent,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          ket_prim,ket_gen,ket_type,ket_cent,
c ...........                             canAB   canCD   canPQ
     &          dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........   DryRun
     &          .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &          bra_rel,ket_rel,ss_one_cent,do_ssss,0)
          else
            call hf2(dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &          dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
     &          .true.,dbl_mb(mb_exndcf(1,1)),memdum)
          end if
*:debug   write(luout,*)' coul ',ibc,ibc,ikc,ikc,memdum
          memsize = max(memsize,memdum)
c
          memdum = memmax
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . exchange
          if (status_rel .and. bra_rel .and. ket_rel) then
            call rel_2e4c_sf (
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          bra_prim,bra_gen,bra_type,bra_cent,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          ket_prim,ket_gen,ket_type,ket_cent,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          bra_prim,bra_gen,bra_type,bra_cent,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          ket_prim,ket_gen,ket_type,ket_cent,
     &          dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
     &          .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &          bra_rel,ket_rel,ss_one_cent,do_ssss,0)
          else
            call hf2(dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &          dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &          dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
     &          .true.,dbl_mb(mb_exndcf(1,1)),memdum)
          end if
*:debug   write(luout,*)' exch ',ibc,ikc,ibc,ikc,memdum
          memsize = max(memsize,memdum)
c
          if (brain.eq.ketin) then
            memdum = memmax
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bb|bk)
            if (status_rel .and. (bra_rel .or. ket_rel)) then
              call rel_2e4c_sf (
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &            bra_rel,bra_rel.and.ket_rel,ss_one_cent,do_ssss,0)
            else
              call hf2(dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum)
            end if
*:debug     write(luout,*)' bbbk ',ibc,ibc,ibc,ikc,memdum
            memsize = max(memsize,memdum)
c
            memdum = memmax
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bk|kk)
            if (status_rel .and. (bra_rel .or. ket_rel)) then
              call rel_2e4c_sf (
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            bra_prim,bra_gen,bra_type,bra_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            ket_prim,ket_gen,ket_type,ket_cent,
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &            bra_rel.and.ket_rel,ket_rel,ss_one_cent,do_ssss,0)
            else
              call hf2(dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &            dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &            dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,
     &            dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
     &            .true.,dbl_mb(mb_exndcf(1,1)),memdum) 
            end if
*:debug     write(luout,*)' bkkk ',ibc,ikc,ikc,ikc,memdum
            memsize = max(memsize,memdum)
*:debug     write(luout,*)' full ',memsize
          endif
00201     continue
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
c
c... add 10% for abbc and abcc class shells not checked for
      memsize = memsize + memsize/10
*      write(luout,'(a,i10)')'emem_2e4c: memsize =',memsize
*      write(luout,'(a,i10)')'emem_2e4c: maxg    =',maxg
c
      end
C>
C> \brief Internal routine to compute the memory required by the 
C> 2-electron 3-center integral code
C>
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_2e3c}
c:tex-This routine computes the memory required for the the two electron
c:tex-three center integrals using the McMurchie-Davidson algorithm. 
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_2e3c(brain, ketin, maxg, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "util.fh"
c::functions
c::passed
c:tex-\begin{verbatim}
      integer brain   !< [Input] bra basis set handle
      integer ketin   !< [Input] ket basis set handle
      integer maxg    !< [Output] max size of 2e integrals buffer
      integer memsize !< [Output] max size of scratch space for 2e integral evaluation
c:tex-\end{verbatim}
c     NB. On input, memsize contains the (dummy) memory limit for 2e integrals
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum, memmax
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
c
#include "bas_exndcf_sfn.fh"
c
      memmax = memsize
      memsize = 0
      maxg = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
*error checking
      if (bra_ucont.eq.0 .or. ket_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' bra handle ',brain
        write(luout,*)' bra name   ',bs_name(bra)
        write(luout,*)' bra contractions ',bra_ucont
        write(luout,*)' ket handle ',ketin
        write(luout,*)' ket name   ',bs_name(ket)
        write(luout,*)' ket contractions ',ket_ucont
        call errquit('emem_2e3c: bad basis info ',911, BASIS_ERR)
      endif
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
          memdum = memmax
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
          call hf2(dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          1, 1, 0, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          bra_prim, bra_gen, bra_type,
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), ket_prim, ket_gen,
     &          ket_type, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          ket_prim, ket_gen, ket_type,
     &          dbl_mb(mb_exndcf(1,1)), 1, .false., .false.,
     &          .false., .true., dbl_mb(mb_exndcf(1,1)),
     &          memdum ) 
          memsize = max(memsize,memdum)
c
00200   continue
00100 continue
c
      call int_nbf_max(brain,max_nbf_bra)
      if (ketin.eq.brain) then
        max_nbf_ket = max_nbf_bra
      else
        call int_nbf_max(ketin,max_nbf_ket)
      endif
      maxg = max_nbf_bra*max_nbf_bra*max_nbf_ket
      maxg = max(maxg,(max_nbf_ket*max_nbf_ket*max_nbf_bra))
c
c... add 10% for safety
      memsize = memsize + memsize/10
*      write(luout,'(a,i10)')'emem_2e3c: memsize =',memsize
*      write(luout,'(a,i10)')'emem_2e3c: maxg    =',maxg
c
      end
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_2e2c}
c:tex-This routine computes the memory required for the the two electron
c:tex-two center integrals using the McMurchie-Davidson algorithm. 
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_2e2c(brain, ketin, maxg, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "util.fh"
c::functions
c::passed
c:tex-\begin{verbatim}
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c:tex-\end{verbatim}
c     NB. On input, memsize contains the (dummy) memory limit for 2e integrals
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum, memmax
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
c
#include "bas_exndcf_sfn.fh"
c
      memmax = memsize
      memsize = 0
      maxg = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
*error checking
      if (bra_ucont.eq.0 .or. ket_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' bra handle ',brain
        write(luout,*)' bra name   ',bs_name(bra)
        write(luout,*)' bra contractions ',bra_ucont
        write(luout,*)' ket handle ',ketin
        write(luout,*)' ket name   ',bs_name(ket)
        write(luout,*)' ket contractions ',ket_ucont
        call errquit('emem_2e2c: bad basis info ',911, BASIS_ERR)
      endif
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
          memdum = memmax
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
          call hf2(
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), bra_prim, bra_gen,
     &          bra_type, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          1, 1, 0, dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          ket_prim, ket_gen, ket_type,
     &          dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &          dbl_mb(mb_exndcf(1,1)), 1, 1, 0,
     &          dbl_mb(mb_exndcf(1,1)), 1, .false., .false.,
     &          .false., .true., dbl_mb(mb_exndcf(1,1)),
     &          memdum)
*          write(luout,*)' emem_2e2c: memdum',memdum, memsize
          memsize = max(memsize,memdum)
c
00200   continue
00100 continue
c
      call int_nbf_max(brain,max_nbf_bra)
      if (ketin.eq.brain) then
        max_nbf_ket = max_nbf_bra
      else
        call int_nbf_max(ketin,max_nbf_ket)
      endif
      maxg = max_nbf_bra*max_nbf_ket
c
c... add 10% for safety
      memsize = memsize + memsize/10
*      write(luout,'(a,i10)')'emem_2e2c: memsize =',memsize
*      write(luout,'(a,i10)')'emem_2e2c: maxg    =',maxg
c
      end
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_2e4c\_full}
c:tex-This routine computes the memory required by the McMurchie-Davidson 
c:tex-algorithm for 4 center two electron integrals based on the 
c:tex-basis sets used.  
c:tex-
c:tex-The exact algorithm is an N$^4$ and costly.  This routine is used
c:tex-primarily as a debugging tool
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_2e4c_full(brain, ketin, maxg, memsize)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
c::functions
      logical cando_sp_u
      external cando_sp_u
c::passed
c:tex-\begin{verbatim}
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c:tex-\end{verbatim}
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, jbc, kbc, lbc
      integer kbc_end
      integer memdum
      integer i_type, j_type, k_type, l_type
      integer i_prim, j_prim, k_prim, l_prim
      integer i_gen,  j_gen,  k_gen,  l_gen
      integer i_cent, j_cent, k_cent, l_cent
      integer max_nbf_bra,max_nbf_ket
c
      integer approx_maxg
      integer approx_memsize
c
      logical bra_cando_sp, bra_rel
      logical ket_cando_sp, ket_rel
      logical status_rel
c
      double precision factor
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      call emem_2e4c(brain,ketin,approx_maxg, approx_memsize)
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
      write(luout,*)' bra basis ',bs_name(bra)(1:len_bs_name(bra))
      write(luout,*)' ket basis ',bs_name(ket)(1:len_bs_name(ket))
c
*error checking
      if (bra_ucont.eq.0 .or. ket_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' bra handle ',brain
        write(luout,*)' bra name   ',bs_name(bra)
        write(luout,*)' bra contractions ',bra_ucont
        write(luout,*)' ket handle ',ketin
        write(luout,*)' ket name   ',bs_name(ket)
        write(luout,*)' ket contractions ',ket_ucont
        call errquit('emem_2e4c_full: bad basis info ',911, BASIS_ERR)
      endif
c 
c relativistic?
c
      status_rel = dyall_mod_dir .and. (.not.nesc_1e_approx)
     &    .and. (brain .eq. ketin) .and. (brain .eq. ao_bsh)
c
      do 00100 ibc = 1, bra_ucont
        i_type = infbs_cont(CONT_TYPE,ibc,bra)
        i_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        i_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        if (i_type.lt.0) then
          i_type = abs(i_type)
          i_gen  = 1
        endif
        i_cent  = (sf_ibs_cn2ce(ibc,bra))
c
        do 00200 jbc = 1, ibc
          j_type = infbs_cont(CONT_TYPE,jbc,bra)
          j_prim = infbs_cont(CONT_NPRIM,jbc,bra)
          j_gen  = infbs_cont(CONT_NGEN,jbc,bra)
          if (j_type.lt.0) then
            j_type = abs(j_type)
            j_gen  = 1
          endif
          j_cent  = (sf_ibs_cn2ce(jbc,bra))
c
          bra_cando_sp = cando_sp_u(brain,ibc,jbc)
          bra_rel = infbs_cont(CONT_RELLS,ibc,bra) .ne. 0
     &        .and. infbs_cont(CONT_RELLS,jbc,bra) .ne. 0
c
          if (bra.eq.ket) then
            kbc_end = jbc
          else
            kbc_end = ket_ucont
          endif
c
          do 00300 kbc = 1, kbc_end
            k_type = infbs_cont(CONT_TYPE,kbc,ket)
            k_prim = infbs_cont(CONT_NPRIM,kbc,ket)
            k_gen  = infbs_cont(CONT_NGEN,kbc,ket)
            if (k_type.lt.0) then
              k_type = abs(k_type)
              k_gen  = 1
            endif
            k_cent  = (sf_ibs_cn2ce(kbc,ket))
c
            do 00400 lbc = 1, kbc
              l_type = infbs_cont(CONT_TYPE,lbc,ket)
              l_prim = infbs_cont(CONT_NPRIM,lbc,ket)
              l_gen  = infbs_cont(CONT_NGEN,lbc,ket)
              if (l_type.lt.0) then
                l_type = abs(l_type)
                l_gen  = 1
              endif
              l_cent  = (sf_ibs_cn2ce(lbc,ket))
c
              ket_cando_sp = cando_sp_u(ketin,kbc,jbc)
              ket_rel = infbs_cont(CONT_RELLS,kbc,ket) .ne. 0
     &            .and. infbs_cont(CONT_RELLS,lbc,ket) .ne. 0
              if (bra_cando_sp.and.ket_cando_sp) goto 00201
c
              memdum = 14 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
              if (status_rel .and. (bra_rel .or. ket_rel)) then
                call rel_2e4c_sf (
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              i_prim, i_gen, i_type, i_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              j_prim, j_gen, j_type, j_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              k_prim, k_gen, k_type, k_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              l_prim, l_gen, l_type, l_cent,
c ...........                                 canAB   canCD   canPQ
     &              dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........       DryRun
     &              .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &              bra_rel,ket_rel,ss_one_cent,do_ssss,0)
              else
                call hf2(
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), i_prim, i_gen, i_type,
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), j_prim, j_gen, j_type,
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), k_prim, k_gen, k_type,
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), l_prim, l_gen, l_type,
     &              dbl_mb(mb_exndcf(1,1)), 1, .false., .false.,
     &              .false.,.true., dbl_mb(mb_exndcf(1,1)),
     &              memdum)
              end if
#if defined(API_DEBUG)
              write(luout,*)' coul ',ibc,jbc,kbc,lbc,memdum
#endif
              memsize = max(memsize,memdum)
c
              memdum = 14 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . exchange
              if (status_rel .and. (bra_rel .or. ket_rel)) then
                call rel_2e4c_sf (
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              i_prim, i_gen, i_type, i_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              l_prim, l_gen, l_type, l_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              k_prim, k_gen, k_type, k_cent,
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &              j_prim, j_gen, j_type, j_cent,
c ...........                                 canAB   canCD   canPQ
     &              dbl_mb(mb_exndcf(1,1)),1,.false.,.false.,.false.,
c ...........       DryRun
     &              .true.,dbl_mb(mb_exndcf(1,1)),memdum,
     &              bra_rel,ket_rel,ss_one_cent,do_ssss,0)
              else
                call hf2(
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), i_prim, i_gen, i_type,
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), l_prim, l_gen, l_type,
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), k_prim, k_gen, k_type,
     &              dbl_mb(mb_exndcf(1,1)), dbl_mb(mb_exndcf(1,1)),
     &              dbl_mb(mb_exndcf(1,1)), j_prim, j_gen, j_type,
     &              dbl_mb(mb_exndcf(1,1)), 1, .false., .false.,
     &              .false., .true., dbl_mb(mb_exndcf(1,1)),
     &              memdum)
              end if
              memsize = max(memsize,memdum)
#if defined(API_DEBUG)
              write(luout,*)' exch ',ibc,jbc,kbc,lbc,memdum
              write(luout,*)' full ',memsize
#endif
00201         continue
00400       continue
00300     continue
00200   continue
00100 continue
c
      call int_nbf_max(brain,max_nbf_bra)
      if (ketin.eq.brain) then
        max_nbf_ket = max_nbf_bra
      else
        call int_nbf_max(ketin,max_nbf_ket)
      endif
      maxg = max_nbf_bra*max_nbf_bra*max_nbf_ket*max_nbf_ket
c
      factor = dble(maxg)/dble(approx_maxg)
      write(luout,'(a,i10)')'emem_2e4c     : maxg    =',approx_maxg
      write(luout,'(a,i10)')'emem_2e4c_full: maxg    =',maxg
      write(luout,'(a,f10.3)')'             maxg factor =',factor
      factor = dble(memsize)/dble(approx_memsize)
      write(luout,'(a,i10)')'emem_2e4c     : memsize =',approx_memsize
      write(luout,'(a,i10)')'emem_2e4c_full: memsize =',memsize
      write(luout,'(a,f10.3)')'          memsize factor =',factor
      call util_flush(6)
c
      end
***********************************************************************
c:tex-% Part of the API internal routines.
c:tex-\subsection{int\_nbf\_max}
c:tex-This routine computes the maximum cartesian nbf for a given basis
c:tex-set.  Used in many memory computing routines to deterine maximum
c:tex-buffer sizes needed for integral computations.  This also includes
c:tex-any general contraction information.
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_nbf_max(basisin,nbf_max)
c:tex-\end{verbatim}
      implicit none
c
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "stdio.fh"
#include "int_nbf.fh"
#include "util.fh"
c::passed
c:tex-\begin{verbatim}
      integer basisin  ! [input] basis set handle
      integer nbf_max  ! [output] maximum number of basis functions
c:tex-\end{verbatim}
c::local
      integer nbf
      integer basis
      integer itype
      integer igen
      integer ucon
      integer num_ucon
c
      basis    = basisin + BASIS_HANDLE_OFFSET
      num_ucon = infbs_head(Head_NCONT,basis)
c
      if (num_ucon.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' basis handle ',basisin, basis
        write(luout,*)' basis name   ',bs_name(basis)
        write(luout,*)' basis contractions ',num_ucon
        call errquit('int_nbf_max: bas basis info ',911, BASIS_ERR)
      endif
c
      nbf_max = 0
      do ucon = 1,num_ucon
        itype = infbs_cont(Cont_TYPE,ucon,basis)
        igen  = infbs_cont(Cont_NGEN,ucon,basis)
        if (itype.lt.0) then
          nbf   = int_nbf_x(itype)
        else
          nbf   = int_nbf_x(itype)*igen
        endif
        nbf_max = max(nbf_max,nbf)
      enddo
      end
c:tex-% Part of the API internal routines.
c:tex-\subsection{int\_mem\_zero}
c:tex-This routine zeros the memory pointers in the apiP.fh common
c:tex-that pertain to the memory utilization of the integral suite.
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_mem_zero()
c:tex-\end{verbatim}
c:tex-There are no formal arguments to this routine
c:tex-
      implicit none
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
      mem_2e4c   = 0   ! scr memory size for 4 center integrals
      isz_2e4c   = 0   ! max buffer size for 4 center integrals
      mem_1e     = 0   ! scr memory size for 1 electron integrals
      isz_1e     = 0   ! max buffer size for 1 electron integrals
      mem_1e_h1  = 0   ! scr memory size for h1-1 electron integrals
      isz_2e3c   = 0   ! max buffer size for 3 center integrals
      mem_2e3c   = 0   ! scr memory size for 3 center integrals
      isz_2e2c   = 0   ! max buffer size for 2 center integrals
      mem_2e2c   = 0   ! scr memory size for 2 center integrals
      isz_3ov    = 0   ! max buffer size for 3 center overlap integrals
      mem_3ov    = 0   ! scr memory size for 3 center overlap integrals
      iszb_2e4c  = 0   ! max buffer size for blocked 4 center integrals
      memb_2e4c  = 0   ! scr memory size for blocked 4 center integrals
      init_int   = 0   ! integral initialization flag
      init_intd  = 0   ! derivative initialization flag
      init_intdd = 0   ! 2nd derivative initialization flag
c
      end
c:tex-% Part of the API internal routines.
c:tex-\subsection{api\_is\_ecp\_basis}
c:tex-This routine identifies basis set handles to INT-API 
c:tex-that store ECP information.
c:tex-This routine should never be called directly by a 
c:tex-NWChem application module.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function api_is_ecp_basis(basisin)
c:tex-\end{verbatim}
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
c:tex-\begin{verbatim}
      integer basisin ! [input] basis set handle
c:tex-\end{verbatim}
c:tex-Return value is true if ``basisin'' represents an ECP
c:tex-
c
      integer basis
c
      basis = basisin + basis_handle_offset
      if (infbs_head(Head_ECP,basis).eq.1) then
        api_is_ecp_basis = .true.
      else
        api_is_ecp_basis = .false.
      endif
      end
      logical function api_is_atomic_basis(basisin)
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin ! [input] basis set handle
* returns true if basis set is an atomic style basis
      integer basis
c
      basis = basisin + basis_handle_offset
      if (infbs_head(Head_ECP,basis).eq.0) then
        api_is_atomic_basis = .true.
      else
        api_is_atomic_basis = .false.
      endif
      end
      logical function api_is_so_basis(basisin)
      implicit none
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
      integer basisin ! [input] basis set handle
* returns true if basis set is an so style basis
      integer basis
c
      basis = basisin + basis_handle_offset
      if (infbs_head(Head_ECP,basis).eq.2) then
        api_is_so_basis = .true.
      else
        api_is_so_basis = .false.
      endif
      end
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{emem\_1e\_pvp}
c:tex-This routine computes the memory for the p.Vp type 1e
c:tex-integrals based on the basis sets used. These integrals are
c:tex-in essence double derivative potential energy integrals of which
c:tex-only the dot product (diagonal) integrals are included.
c:tex-This routine should never be called directly by a
c:tex-NWChem application module.
c:tex-
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine emem_1e_pvp(brain, ketin, max1e, memsize, rel_typ)
c:tex-\end{verbatim}
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
#include "util.fh"
c::passed
c:tex-\begin{verbatim}
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer rel_typ ! [input] type of integrals to be computed
      integer max1e   ! [output] max size of 1e integrals buffer
      integer memsize ! [output] max size of scratch space for 1e integral evaluation
c:tex-\end{verbatim}
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
      integer rel_dbg
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
c     rel_typ = 1
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
*error checking
      if (bra_ucont.eq.0 .or. ket_ucont.eq.0) then
        write(luout,*)' the number of unique contractions is zero'
        write(luout,*)' bra handle ',brain
        write(luout,*)' bra name   ',bs_name(bra)
        write(luout,*)' bra contractions ',bra_ucont
        write(luout,*)' ket handle ',ketin
        write(luout,*)' ket name   ',bs_name(ket)
        write(luout,*)' ket contractions ',ket_ucont
        call errquit('emem_1e_rel: bad basis info ',911, BASIS_ERR)
      endif
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
          memdum = 14 000 000
*
          call rel_pvp(dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &         dbl_mb(mb_exndcf(1,1)),bra_prim,bra_gen,bra_type,1,
     &         dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &         dbl_mb(mb_exndcf(1,1)),ket_prim,ket_gen,ket_type,1,
     &         dbl_mb(mb_exndcf(1,1)),dbl_mb(mb_exndcf(1,1)),
     &         dbl_mb(mb_exndcf(1,1)),1,
     &         dbl_mb(mb_exndcf(1,1)),1,.false.,.true.,
     &         dbl_mb(mb_exndcf(1,1)),memdum,'emem_1e_pvp',
     &         rel_dbg,rel_typ)
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
      max1e = max_nbf_bra*max_nbf_ket
c... insure scratch is large enough to hold 3*max1e
      memsize = max(3*max1e,memsize)
c
#if defined(API_DEBUG)
      write(luout,'(a,i10)')'emem_1e_pvp: memsize =',memsize
      write(luout,'(a,i10)')'emem_1e_pvp: max1e   =',max1e
#endif
c
      end
C> @}
