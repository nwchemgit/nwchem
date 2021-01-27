      subroutine gw_memory(pars)
c
c     Based on tddft_memory
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
   
      double precision ga_dbl, ma_dbl, mo_eris, ao_eris
c
      double precision gw_n3c
      external gw_n3c
c
      ga_dbl = dble(ga_memory_avail())/dble(ma_sizeof(mt_dbl,1,mt_byte))
      ma_dbl = dble(ma_inquire_avail(mt_dbl))
      call ga_dgop(msg_tdmem,ga_dbl,1,'+')
      if (pars%me.eq.0) then
        write(luout,*)
        call util_print_centered(luout,'Memory Information',20,.true.)
        if (ga_uses_ma()) then
          write(luout,9000) int(ga_dbl)
        else
          write(luout,9010) int(ga_dbl),int(ma_dbl)
        endif
      endif
c
 9000 format(10x,'Available GA+MA space size is ',i10,' doubles')
 9010 format(10x,'Available GA space size is    ',i10,' doubles',/,
     $       10x,'Available MA space size is    ',i10,' doubles')
c
c
c
c      ao_eris = gw_n3c(pars)
c
c     3 center ERIs 
c
      mo_eris = dble(nbf_cd * pars%nmo *
     $       max(pars%nocc(1),pars%noqp(1)+pars%nvqp(1)))
      if (pars%me.eq.0) write(luout,9020) int(mo_eris)
 9020 format(10x,'Size of 3-center MO ERI tensor is',i10,' doubles')

      return
      end subroutine
