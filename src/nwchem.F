      program nwchem
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "rtdb.fh"
      integer rtdb
      integer stack
      integer heap
c
c     Create parallel processes and initialize IPC layer
c
      call pbeginf()
c
c     Look for memory directive in the input ... must eventually
c     open the rtdb first so that can get memory directive out of that
c     if it is not in the input
c
      call input_mem_size(stack, heap)
c
c     Initialize local memory allocator
c
      if (.not. ma_init(mt_dbl, stack, heap))
     $     call errquit('nwchem: ma_init failed', -1)
c
c     Intialize global memory allocator
c
      call ga_initialize()
c
c     Read the input and open the RTDB
c
      call start(rtdb)
c
c     Call driver routine to do something useful
c
      call control(rtdb)
c
c     Close the RTDB
c
      if (.not. rtdb_par_close(rtdb, 'keep'))
     $     call errquit('nwchem: rtdb_close failed', rtdb)
c
c     Tidy up the global memory
c
      call ga_terminate()
c
c     Tidy up the parallel environment
c
      call pend()
c
c     Done
c
      end
      subroutine start(rtdb)
      implicit none
#include "mafdecls.h"
#include "rtdb.fh"
      integer rtdb              ! [output]
      character*255 rtdb_name
      character*255 file_prefix
      logical startup
c
c     Figure out if this is a start/restart and what the database
c     and file prefix are
c
      call input_start_opt(rtdb_name, file_prefix, startup)
c
c     Actually open the database and insert the file_prefix
c
      if (startup) then
         if (.not. rtdb_par_open(rtdb_name, 'empty', rtdb))
     $        call errquit('start: rtdb_open empty failed ', 0)
      else
         if (.not. rtdb_par_open(rtdb_name, 'old', rtdb))
     $        call errquit('start: rtdb_open old failed ', 0)
      endif
      if (.not. rtdb_cput(rtdb, 'file_prefix', 1, file_prefix))
     $     call errquit('start: rtdb_put of file_prefix failed', 0)
c      
c     Make sure input is all ready to go
c
      rewind 5
      call inp_init(5,6)
c
c     Parse input data and shove into the database
c
      call input_parse(rtdb)
c
c     Done
c
      end
      subroutine control(rtdb)
      implicit none
      integer rtdb              ! [input] data base handle
#include "rtdb.fh"
c
      call print_integs(rtdb, 'geometry','mo basis')
c
      if (.not. rtdb_print(rtdb, .true.))
     $     call errquit('control: rtdb_print failed', 0)
c
      end
      subroutine print_integs(rtdb, geom_name, basis_name)
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "inp.fh"
      integer rtdb
      character*(*) geom_name, basis_name
      integer geom
      integer basis
      integer ncenter
      integer icent, jcent, kcent, lcent
      integer ish, jsh, ksh, lsh
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer i, j, k, l, ijkl
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(256), xi(3), xj(3), xk(3), xl(3)
      double precision qi, qj, qk, ql
      character*16 tagi, tagj, tagk, tagl
c
      if (.not. geom_create(geom, geom_name))
     $     call errquit('print_integs: geom_create failed?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
     $     call errquit('print_integs: geom_load failed', 0)
      if (.not. geom_print(geom))
     $     call errquit('print_integs: geom_print failed', 0)
      if (.not. bas_create(basis, basis_name))
     $     call errquit('print_integs: basis create failed', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
     $      call errquit('print_integs: basis load failed', 0)
      if (.not. bas_print(basis))
     $      call errquit('print_integs: basis print failed', 0)
      if (.not. gbs_map_print(basis))
     $     call errquit('print_integs: gbs_map_print failed', 0)
c
      write(6,2) geom_name(1:inp_strlen(geom_name)),
     $     basis_name(1:inp_strlen(basis_name))
 2    format(/' Integrals for geometry "',a,'" with basis "',a,'"'/)
c
c     Partially redundant four-fold loop thru centers
c
      if (.not. geom_ncent(geom, ncenter))
     $     call errquit('print_integs: geom_ncent failed', 0)
      do icent = 1, ncenter
         do jcent = 1, icent
            do kcent = 1, icent
               do lcent = 1, kcent
c
c     Figure out coords and contractions of each center
c
                  if (.not. geom_cent_get(geom, icent, tagi, xi, qi))
     $                 call errquit('print_integs: cent_get ?', icent)
                  if (.not. geom_cent_get(geom, jcent, tagj, xj, qj))
     $                 call errquit('print_integs: cent_get ?', jcent)
                  if (.not. geom_cent_get(geom, kcent, tagk, xk, qk))
     $                 call errquit('print_integs: cent_get ?', kcent)
                  if (.not. geom_cent_get(geom, lcent, tagl, xl, ql))
     $                 call errquit('print_integs: cent_get ?', lcent)
c                  
                  if (.not. bas_ce2cnr(basis, icent, ishlo, ishhi))
     $                 call errquit('print_integs: ce2cnr ?', icent)
                  if (.not. bas_ce2cnr(basis, jcent, jshlo, jshhi))
     $                 call errquit('print_integs: ce2cnr ?', jcent)
                  if (.not. bas_ce2cnr(basis, kcent, kshlo, kshhi))
     $                 call errquit('print_integs: ce2cnr ?', kcent)
                  if (.not. bas_ce2cnr(basis, lcent, lshlo, lshhi))
     $                 call errquit('print_integs: ce2cnr ?', lcent)
c
c     Partially redundant four-fold loop thru contractions
c
                  do ish = ishlo, ishhi
                     do jsh = jshlo, jshhi
                        do ksh = kshlo, kshhi
                           do lsh = lshlo, lshhi
c
c     Now compute the SP integrals over these shells
c
                              call genr70(
     $                             basis, ish, xi, jsh, xj,
     $                             basis, ksh, xk, lsh, xl,
     $                             eri)
c
c     Print the integrals out
c 
                              if (.not. bas_cn2bfr(basis, ish, ilo,ihi))
     $                             call errquit('print_integ: cn2bf?',0)
                              if (.not. bas_cn2bfr(basis, jsh, jlo,jhi))
     $                             call errquit('print_integ: cn2bf?',0)
                              if (.not. bas_cn2bfr(basis, ksh, klo,khi))
     $                             call errquit('print_integ: cn2bf?',0)
                              if (.not. bas_cn2bfr(basis, lsh, llo,lhi))
     $                             call errquit('print_integ: cn2bf?',0)
c
                              ijkl = 0
                              do i = 1, ilo, ihi
                                 do j = jlo, jhi
                                    do k = klo, khi
                                       do l = llo, lhi
                                          ijkl = ijkl + 1
                                          write(6,1) i,j,k,l,eri(ijkl)
                                       enddo
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
 1    format(4i6,f20.10)


      if (.not. geom_destroy(geom))
     $     call errquit('print_integs: geom_destroy failed', 0)
      if (.not.bas_destroy(basis))
     $     call errquit('print_integs: basis destroy failed', 0)
c
      end
