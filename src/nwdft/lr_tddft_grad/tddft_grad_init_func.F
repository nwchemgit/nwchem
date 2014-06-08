      subroutine tddft_grad_init_func(numfunc,xfac,cfac,eps)
c
c     $Id$
c
c     The DFT library is a straightforward Fortran90 library without
c     any knowledge of parallelism. In NWChem the initialisation of 
c     the density functional, however, is handled as follows:
c
c     - Processor 0:
c
c       Reads the input and initialises the appropriate data structures
c       When all input has been processed the data structures are 
c       written out to the run-time-data-base (RTDB).
c
c     - All processors:
c
c       Upon entering a Density Functional Theory module read the
c       definition of the functional from the run-time-data-base.
c
c     So the RTDB doubles as a mechanism for communicating the
c     functional to all processors. As a consequence we cannot just
c     modify the input parser to initialise the functional in the 
c     library. Instead a routine is needed that can initialise the
c     library from the data read off the RTDB, which is the task of
c     this routine. 
c
c     Obviously, it is desirable that this routine eventually disappears
c     as it is likely to cause maintenance problems. However, while
c     two functional libraries are being used concurrently it is a 
c     necessary evil.
c
c     Written by Huub van Dam, Sep 2009.
c
      use dftr
      implicit none
#include "stdio.fh"
#include "errquit.fh"
c
c     Input:
c
      integer numfunc                ! the number of different known
                                     ! terms
      double precision xfac(numfunc) ! weights of exchange terms
      double precision cfac(numfunc) ! weights of correlation terms
      double precision eps           ! tolerance on weights
c
c     Local:
c
      integer i
      integer*4 ierror4
      integer nfunc
      parameter (nfunc = 40)
      character*16 xfnc_alias(nfunc), cfnc_alias(nfunc)
      data xfnc_alias / 
     +     'hf_x','s_x','b88_x','hcth_x','b97_x',
     +     'b97_1_x','pbe96_x','gill96_x','pw91_x','hcth120_x',
c 10
     +     'hcth147_x','b98_x','b97g','hcth407_x','hp14',
     +     'opt_x','b97_2_x','pkzb99_x','mpw91_x','407p_x',
c 20
     +     'tpss03_x','b97_3_x','','','ft97_x',
     +     '','','','','rpbe96_x',
c 30
     +     'revpbe96_x','','','','',
     +     '','','','','' /
      data cfnc_alias /
     +     'vwn5_c','lyp88_c','p81_c','p86_c','p91_c',
     +     '','vwn1_rpa_c','vwn1_c','vwn2_c','vwn3_c',
c 10
     +     'vwn4_c','pw91lda_c','hcth_c','b97_c','b97_1_c',
     +     'hcth120_c','hcth147_c','b98_c','b97g_c','hcth407_c',
c 20
     +     'hp14_c','opt_c','b97_2_c','ft97_c','pkzb99_c',
     +     '407p_c','tpss03_c','b97_3_c','','',
c 30
     +     '','','','','',
     +     'be88_c','pb96_c','','','' /
c
c     Code:
c
      if (numfunc.ne.nfunc)
     +  call errquit('tddft_grad_init_func: failed numfunc.ne.nfunc',
     +               numfunc,INPUT_ERR)
c
      call dftr_reset
c
      do i = 1, numfunc
        if (abs(xfac(i)).gt.eps) then
          call dftr_add(xfnc_alias(i),xfac(i),ierror4)
          if (ierror4.eq.dftr_notfound) then
            write(LuOut,*)'*** exchange functional term not found: ',
     +                    xfnc_alias(i)
          endif
        endif
        if (abs(cfac(i)).gt.eps) then
          call dftr_add(cfnc_alias(i),cfac(i),ierror4)
          if (ierror4.eq.dftr_notfound) then
            write(LuOut,*)'*** correlation functional term not found: ',
     +                    cfnc_alias(i)
          endif
        endif
      enddo
c
      end
