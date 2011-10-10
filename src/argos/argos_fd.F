      subroutine argos_fd(isg,xs,fst,fs,iwg,xw,fwt,fw)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer isg(msa),iwg(mwm)
      real*8 xs(msa,3),fst(msa,3),fs(msa,3)
      real*8 xw(mwm,3,mwa),fwt(mwm,3,mwa),fw(mwm,3,mwa)
c
      integer i,j,k,lsa
      real*8 etott,xsk(3),ft(3,3),dx,xdv
c
      if(me.eq.0) write(lfnout,1000)
 1000 format(//,' FINITE DIFFERENCE SOLUTE FORCES',//,
     + '  Atom',
     + '                 Analytic forces       ',
     + '             Finite difference forces  ',
     + '             Deviation                 ',/,
     + '      ',
     + '             fx          fy          fz',
     + '             fx          fy          fz',
     + '            dfx         dfy         dfz',/)
c
      lpair=.false.
      lload=.false.
c      xdv=0.000001
      xdv=dx0sd
c
      do 1 i=1,nsaloc
      fst(i,1)=fs(i,1)
      fst(i,2)=fs(i,2)
      fst(i,3)=fs(i,3)
    1 continue
      do 2 i=1,nwmloc
      do 3 j=1,nwa
      fwt(i,1,j)=fw(i,1,j)
      fwt(i,2,j)=fw(i,2,j)
      fwt(i,3,j)=fw(i,3,j)
    3 continue
    2 continue
      etott=etot
c
      do 4 i=1,nsa
      lsa=0
      ft(1,1)=zero
      ft(2,1)=zero
      ft(3,1)=zero
      ft(1,2)=zero
      ft(2,2)=zero
      ft(3,2)=zero
      ft(1,3)=zero
      ft(2,3)=zero
      ft(3,3)=zero
      do 5 j=1,nsaloc
      if(isg(j).eq.i) then
      lsa=j
      xsk(1)=xs(j,1)
      xsk(2)=xs(j,2)
      xsk(3)=xs(j,3)
      ft(1,1)=fst(j,1)
      ft(2,1)=fst(j,2)
      ft(3,1)=fst(j,3)
      goto 6
      endif
    5 continue
    6 continue
c
      do 7 j=2,3
      do 8 k=1,3
c
      dx=-xdv
      if(j.eq.3) dx=xdv
c
      if(lsa.gt.0) xs(lsa,k)=xs(lsa,k)+dx
c
c     atomic forces and potential energies
c
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call argos_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call argos_prop_proper(0,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
c
      if(lsa.gt.0) then
      xs(lsa,k)=xsk(k)
      ft(k,j)=etott-etot
      endif
c
    8 continue
    7 continue
c
      if(np.gt.1) call ga_dgop(mrg_d14,ft,6,'+')
c
      if(me.eq.0) write(lfnout,1001) i,(ft(j,1),j=1,3),
     + ((ft(j,3)-ft(j,2))/(two*xdv),j=1,3),
     + (ft(j,1)-(ft(j,3)-ft(j,2))/(two*xdv),j=1,3)
 1001 format(i7,5x,3f12.3,3x,3f12.3,3x,3E12.3)
c
    4 continue
c
      return
      end
c $Id$
