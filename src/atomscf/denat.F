      subroutine guess_dens(geom, basis, g_dens)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "cguess.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer geom, basis       ! Handles
      integer g_dens            ! GA returns superposition of atomic densities
c
c..   get starting vectors from superposition of atomic densities
c..
c..   routines involved are (all present in this order) :
c..   denat  : general control routine
c..   datoms : translate orbital info, call atomscf , gather d-matrix
c..   atomd,tracd,trafsd,cmergd,oeigd,teigd,densid,denmad, ..
c..   .. hamild,outpud,shalfd,atcond,starcd,creded,pdfded, ..
c..   .. jacod,orderd,tramad
c..   oeigd has been changed to pick up the 1-elec ints
c..   atcond has been changed to allow for effective nuclear charge
c..   creded/pdfded have been adapted to yield directly a d-matrix
crz   xpsnld: implemented for use of non local pseudo potentials
c..
c...  start all types after one cycle closed shell in routine denscf
c
      integer iprint
      integer iwrite
      logical oprinv
      integer l2
      integer nb, no, ntr, nsq
      integer i10, ipcap, iqcap, ifc, ifo, is, iu, it
      integer ih, idc, idos, idt, idold, iss, ic, icopn, ismin
      integer iqmin, itransf, icc, last, memq, l_scr, k_scr
      data iwrite/6/, iprint/2/
c..
c..    core partitioning
c..
      oprinv = iprint.eq.2
c
c     allow for maximum of 100 bfns on any given atom
c
      ntr = nb*(nb+1)/2
      l2 = ntr
      nsq = nb * nb
c
c     dynamically allocate memory
c
      call guess_mem(memq)
      if (.not. ma_push_get(mt_dbl, memq,'atomic guess', l_scr, k_scr))
     $     call errquit('guess ao dens: ma_push failed', memq)
c
      i10 = k_scr
      ipcap = i10 + l2
      iqcap = ipcap + ntr
      ifc = iqcap + ntr
      ifo = ifc + ntr
      is = ifo + ntr
      iu = is + ntr
      it = iu + ntr
      ih = it + ntr
      idc = ih + ntr
      idos = idc + ntr
      idt = idos + ntr
      idold = idt + ntr
      iss = idold + ntr
      ic = iss + ntr
c
      icopn = ic + nsq
      ismin = icopn + nsq
      iqmin = ismin + nb * no
      itransf = iqmin + nb * no
      icc = itransf + nsq
      last = icc + nsq
      if ((last - i10) .gt. memq) call errquit
     $     ('guess_dens: memory?', memq)
c..
c..   now loop over the atoms / do atomic scf and gather d-matrix
c..
      call datoms(geom, basis, dbl_mb(i10),g_dens,oprinv,
     +     dbl_mb(ipcap), dbl_mb(iqcap), dbl_mb(ifc), dbl_mb(ifo), 
     $     dbl_mb(is), dbl_mb(iu), 
     +     dbl_mb(it), dbl_mb(ih), dbl_mb(idc), dbl_mb(idos), 
     $     dbl_mb(idt), dbl_mb(idold), dbl_mb(iss) ,
     +     dbl_mb(ic), dbl_mb(icopn), dbl_mb(ismin), dbl_mb(iqmin), 
     $     dbl_mb(itransf), dbl_mb(icc), nb) 
c
c..   print if requested
c
      if (oprinv) then
        if (ga_nodeid().eq.0) write (iwrite,6020)
        call ga_print(g_dens)
      endif
c
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('guess ao dens: pop_stack failed', 0)
c
      return
 6020 format (//30x,28('-'),/,30x,'   initial guess density   ',/,30x,
     +        28('-')//)
      end
