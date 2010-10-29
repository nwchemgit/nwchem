      subroutine atscf_atomic_scf(geom, basis, iat, oprint,
     $     nl_out, ncsh_out, nosh_out, nccup_out, nbc_out,
     $     cc_out, energy_out, eps_out,spec)
      implicit none
#include "errquit.fh"
* $Id$
#include "geom.fh"
#include "bas.fh"
#include "cguess.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "stdio.fh"
      integer geom, basis       ! [input] Handles
      integer iat               ! [input] Atom number in NWChem numbering
      logical oprint,spec       ! [input] If true the SCF prints stuff
      integer nl_out            ! [output] No. of angular momenta (= sym)
      integer ncsh_out(*)       ! [output] No. of closed shells per sym
      integer nosh_out(*)       ! [output] No. of open shells per sym
      integer nccup_out(*)      ! [output] Occupation no. of open shells
      integer nbc_out(*)        ! [output] No. of basis shells per sym
      double precision cc_out(*) ! [output] MO coefficients
      double precision energy_out ! [output] Energy
      double precision eps_out(*) ! [output] Orbital energies
c
c     Perform an atomic SCF for a given atom and return energetics
c     and enough information to reconstruct the wavefunction
c
      integer l2
      integer ntr, nsq
      integer i10, ipcap, iqcap, ifc, ifo, is, iu, it, i
      integer ih, idc, idos, idt, idold, iss, ic, icopn, ismin
      integer iqmin, itransf, icc, last, memq, l_scr, k_scr
c
      integer ii, ii_out, ieps, ieps_out
c
c     allow for maximum of 100 bfns (NWChem shells) on any given atom
c
      ntr = nb*(nb+1)/2
      l2 = ntr
      nsq = nb * nb
c
c     dynamically allocate memory
c
      call guess_mem(memq)
      if (.not. ma_push_get(mt_dbl, memq,'atomic guess', l_scr, k_scr))
     $     call errquit('atomscf:guess_dens: ma_push failed', memq,
     &       MA_ERR)
      call dfill(memq,0.0d00,dbl_mb(k_scr),1)
c
c     This is verbatim from denat
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
     $     ('atomscf:guess_dens: memory?', memq, MEM_ERR)
c..
      call atscf_do_atscf(geom, basis, iat,
     $     dbl_mb(i10), oprint,
     +     dbl_mb(ipcap), dbl_mb(iqcap), dbl_mb(ifc), dbl_mb(ifo), 
     $     dbl_mb(is), dbl_mb(iu), 
     +     dbl_mb(it), dbl_mb(ih), dbl_mb(idc), dbl_mb(idos), 
     $     dbl_mb(idt), dbl_mb(idold), dbl_mb(iss) ,
     +     dbl_mb(ic), dbl_mb(icopn), dbl_mb(ismin), dbl_mb(iqmin), 
     $     dbl_mb(itransf), dbl_mb(icc), nb,spec)
c
c     Copy stuff into output arguments
c
      energy_out = energ
      nl_out = nsym             ! Highest basis ang. mom. (s=1,p=2,d=3,f=4)
      ii     = icc
      ii_out = 1
      ieps     = 1
      ieps_out = 1
      do i = 1, nsym
         ncsh_out(i) = ncsh(i)
         nosh_out(i) = nosh(i)
         nccup_out(i)= nccup(i)
         nbc_out(i)  = nbc(i)
c     Copy only occupied orbitals and eigenvalues
         call dcopy(nbc(i)*(ncsh(i)+nosh(i)),dbl_mb(ii),1,
     $        cc_out(ii_out),1)
         ii     = ii     + nbc(i)*nbc(i)
         ii_out = ii_out + nbc(i)*(ncsh(i)+nosh(i))
         call dcopy(ncsh(i)+nosh(i),eps(ieps),1,eps_out(ieps_out),1)
         ieps     = ieps     + ncsh(i)+nosh(i)
         ieps_out = ieps_out + ncsh(i)+nosh(i)
      enddo 
c
      if (oprint) then
c
c     This only to debug this interface
c
         write(6,*) ' Atomic SCF for atom ', iat
         write(6,*) ' nl_out    ', nl_out
         write(6,*) ' ncsh   ', (ncsh_out(i),i=1,nl_out)
         write(6,*) ' nosh   ', (nosh_out(i),i=1,nl_out)
         write(6,*) ' nccup  ', (nccup_out(i),i=1,nl_out)
         write(6,*) ' nbc    ', (nbc_out(i),i=1,nl_out)
         write(6,*) ' energy ', energy_out
         ii_out   = 1
         ieps_out = 1
         do i = 1, nl_out
            write(6,*) ' Occupied eigenvalues for symmetry ', i
            call output(eps_out(ieps_out),1,(ncsh_out(i)+nosh_out(i)),
     $           1,1,nbc_out(i),1,1)
            ieps_out = ieps_out + (ncsh_out(i)+nosh_out(i))
            write(6,*) ' Occupied eigenvectors for symmetry ', i
            call output(cc_out(ii_out),1,nbc_out(i),
     $           1,(ncsh_out(i)+nosh_out(i)), nbc_out(i),nbc_out(i),1)
            ii_out = ii_out + nbc(i)*(ncsh(i)+nosh(i))
         enddo
      endif
c
      if (.not. ma_pop_stack(l_scr))
     $     call errquit
     &    ('atomscf:guess_dens: pop_stack failed', 0, MEM_ERR)
c
      end
      subroutine atscf_do_atscf(geom, basis, iat,
     $     hatom, oprin,
     +     pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     +     cvec, copn, smin, qmin, transf, cc, nbb,spec)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "cguess.fh"
#include "msgids.fh"
      integer geom, basis
      integer iat
c     
      integer nbb               ! No. of primtives used to dimension arrays
      real *8 hatom(*)
      real *8 pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      real *8 h(*), dc(*), dos(*), dt(*), dold(*), ss(*)
      real *8 cvec(*), copn(*), smin(nbb,*), qmin(nbb,*),transf(*),cc(*)
      logical oprin,spec
      integer iwrite
      character*16 tag, element
c
      logical ospherical
c     
      integer ic(4,nb),iiloc(nb,4)
c     
      integer isymax
      real *8 znps
      data iwrite/6/
c     
      call atscf_bas_gather(geom, basis, iat, ic, iiloc,
     $     isymax, hatom, znps, element, tag, ospherical)
c     
      call atomd(oprin,iwrite,znps,ic,isymax,hatom,
     +     pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     +     cvec, copn, smin, qmin, transf, cc, nbb,
     +     element, tag,spec)
c
      end
