      subroutine intd_init(rtdb,nbas,bases)
c $Id: intd_init.F,v 1.6 1996-07-11 22:21:19 d3e129 Exp $
c
c Initializes integral derivative code 
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "apiP.fh"
#include "rtdb.fh"
c::passed
      integer rtdb        ! run time data base handle
      integer nbas        ! number of basis sets to be used
      integer bases(nbas) ! basis set handles
c::local
      integer ibas, ang2use, angm
      logical status
      integer intd_memtmp
c
*      call int_init(nbas,bases)
c
      if(init_intd.eq.1) then
        write(6,*)' warning nested intd_inits'
        write(6,*)' intd_init already called '
        call util_flush(6)
      endif
      status = .true.
      do 00100 ibas=1,nbas
        status = status .and. bas_check_handle(bases(ibas),'intd_init')
00100 continue

      if (.not.status) then
        write(6,*)' at least one basis handle not valid'
        do 00200 ibas = 1,nbas
          write(6,'(a,i5)')
     &           ' basis set handle ',bases(ibas)
00200   continue
        call errquit('intd_init: basis handles hosed ',nbas)
      endif
*      write(6,*)' intd_init: basis set handles valid '
c
c initialize defnxyz routines
c      
      ang2use = -1
      do 00300 ibas = 1,nbas
        if(.not.bas_high_angular(bases(ibas),angm))
     &         call errquit('intd_init: angm error',angm)
        ang2use = max(ang2use,angm)
00300 continue
c.. for first derivatives add 1
      call defNxyz(ang2use+1)
c.. read in approximate memory from rtdb
c.  
c.. parameter is default value used in hf2d
c.. input memory in words should be scaled by 1/12 for same range
      intd_memthresh = intd_memp
      if (rtdb_get(rtdb,'intd:approxmem',MT_INT,1,intd_memtmp)) then
        if(ga_nodeid().eq.0)then
          write(6,'(/a,i10,a/)')
     &        ' approximate memory for derivative integrals set to:',
     &        intd_memtmp
        endif
        intd_memthresh = intd_memtmp/12
      endif
c
c... generate memory requirements and store in structures in apiP.fh
c
      call exactd_mem(bases,nbas)
      call sp_der_init(nbas,bases)
      call int_acc_std()
      init_int = 1
      init_intd = 1
      end
      subroutine sp_der_init(nbas,bases)
      implicit none
#include "apispP.fh"
#include "bas.fh"
#include "mafdecls.fh"
c::passed
      integer nbas
      integer bases(nbas)
c::local
      integer ncont,icont,ibas
      integer mytype,nprim,ngen,sphorcart
c
      no_sp_conts = .true.
      do 00100 ibas = 1,nbas
        if (.not.bas_numcont(bases(ibas),ncont))
     &         call errquit('sp_der_init: bas_numcont failed',911)
        do 00200 icont = 1,ncont
          if(.not.bas_continfo
     &           (bases(ibas),icont,mytype,nprim,ngen,sphorcart))
     &           call errquit('sp_der_init: bas_continfo failed',911)

          if (mytype.lt.0) then
            no_sp_conts = .false.
            goto 00201
          endif
00200   continue
00100 continue
00201 continue
c
      if (no_sp_conts) return
c
      write(6,*)' sp functions cannot be used in the integral ',
     &       'derivative code yet'
      call errquit('sp_der_init',911)
      end

      
