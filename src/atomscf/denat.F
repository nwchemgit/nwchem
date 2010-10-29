      subroutine guess_dens(rtdb, geom, basis, g_dens)
      implicit none
      integer rtdb, geom, basis, g_dens
c
      logical odens
      double precision eatom
c
      odens = .true.
c
      call guess_dens_gen(rtdb, geom, basis, odens, g_dens, eatom)
c
      end
      subroutine guess_dens_gen(rtdb, geom, basis, odens, g_dens, eatom)
C$Id$
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "cguess.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "atomz.fh"
#include "inp.fh"
#include "stdio.fh"
      integer rtdb              ! Run Time Data Base Handle
      integer geom, basis       ! Handles
      logical odens             ! [input] If true make the density
      integer g_dens            ! GA returns superposition of atomic densities
      double precision eatom    ! returns sum of atomic SCF energies
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
      integer iwrite
      logical oprint_dens, oprint_scf, oprint_ener, oprint_mod_z,spec
      integer l2, iii, iiilen
      integer ntr, nsq
      integer i10, ipcap, iqcap, ifc, ifo, is, iu, it
      integer ih, idc, idos, idt, idold, iss, ic, icopn, ismin
      integer iqmin, itransf, icc, last, memq, l_scr, k_scr
c
      integer ma_type
      character*26 date
c
      data iwrite/6/, oprint_dens/.false./, oprint_scf/.false./
c
      oprint_dens = (odens .and. 
     $     util_print('atomic guess density',print_debug))
      oprint_scf  = util_print('atomic scf', print_debug)
      oprint_ener = util_print('mo guess', print_default)
      oprint_mod_z = util_print('mo guess', print_default)
c
      iwrite = LuOut
      if (ga_nodeid().eq.0 .and. 
     $     (oprint_dens .or. oprint_scf .or. oprint_ener)) then
         write(LuOut,*)
         call util_print_centered(LuOut,
     $        'Superposition of Atomic Density Guess', 23, .true.)
         write(LuOut,*)
      endif
c..
c..    core partitioning
c..
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
     $     call errquit('guess_dens: ma_push failed', memq, MA_ERR)
      call dfill(memq,0.0d00,dbl_mb(k_scr),1)
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
     $     ('guess_dens: memory?', memq, MEM_ERR)
c..
c.. parse data from rtdb about atoms
c..
*      if (.not.rtdb_get(rtdb,'atomscf:num_z',mt_int,1,num_z))
*     &    num_z = 0
c
      if (.not. rtdb_get(rtdb, 'atomscf:special', mt_log, 1, spec))
     $     spec = .false.
c
      if (.not. rtdb_get_info(rtdb, 'atomscf:z', ma_type, num_z, date))
     $     num_z = 0
c
      if (num_z .gt. max_num_z)
     &    call errquit
     &    ('guess_dens: max_num_z < num_z  delta =',(max_num_z-num_z),
     &       UNKNOWN_ERR)
      if (num_z .gt. 0) then
        if (.not.rtdb_get(rtdb,'atomscf:z',mt_dbl,num_z,z))
     &      call errquit('guess_dens: read of atomscf:z failed ',911,
     &       UNKNOWN_ERR)
        if (.not.rtdb_cget(rtdb,'atomscf:tags_z',num_z,tags_z))
     &      call errquit
     &      ('guess_dens: read of atomscf:tag_z failed',911,
     &       INPUT_ERR)
        if (oprint_mod_z.and.(ga_nodeid().eq.0)) then
          do iii = 1,num_z
            iiilen = inp_strlen(tags_z(iii))
            write(luout,'(1x,a,f7.2,a,a)')
     &          'Atomic density modified by a delta charge of',
     &          z(iii),' for atom tag ',
     &          tags_z(iii)(1:iiilen)
          enddo
          write(luout,*)' '
        endif
      endif
c..
c..   now loop over the atoms / do atomic scf and gather d-matrix
c..
      call datoms(geom, basis, dbl_mb(i10),g_dens,oprint_scf,
     +     dbl_mb(ipcap), dbl_mb(iqcap), dbl_mb(ifc), dbl_mb(ifo), 
     $     dbl_mb(is), dbl_mb(iu), 
     +     dbl_mb(it), dbl_mb(ih), dbl_mb(idc), dbl_mb(idos), 
     $     dbl_mb(idt), dbl_mb(idold), dbl_mb(iss) ,
     +     dbl_mb(ic), dbl_mb(icopn), dbl_mb(ismin), dbl_mb(iqmin), 
     $     dbl_mb(itransf), dbl_mb(icc), nb, oprint_ener,
     $     odens, eatom,spec) 
c
c..   print if requested
c
      if (oprint_dens) then
        if (ga_nodeid().eq.0) write (iwrite,6020)
        call ga_print(g_dens)
      endif
c
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('guess_dens: pop_stack failed', 0, MEM_ERR)
c     
      if (util_print('debugatscf',print_never))
     $     call atscf_debug(geom,basis,spec)
c
      return
 6020 format (//30x,28('-'),/,30x,'   initial guess density   ',/,30x,
     +        28('-')//)
      end
      subroutine atscf_debug(geom,basis,spec)
      implicit none
#include "errquit.fh"
      integer geom,basis
      logical spec
#include "geom.fh"
c
      integer nl_out, iat, nat
      integer ncsh_out(4),nosh_out(4),nccup_out(4),nbc_out(4)
      double precision cc_out(2000), energy_out, eps_out(50)
c
      if (.not. geom_ncent(geom, nat)) call errquit('atscf_debug',0,
     &       UNKNOWN_ERR)
      do iat = 1, nat
         call atscf_atomic_scf(geom, basis, iat, .true.,
     $        nl_out, ncsh_out, nosh_out, nccup_out, nbc_out,
     $        cc_out, energy_out, eps_out,spec)
      enddo
c
      end

