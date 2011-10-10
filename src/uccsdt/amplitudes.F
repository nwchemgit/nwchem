      logical function uccsdt_ampfile_init(rtdb)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "dra.fh"
#include "rtdb.fh"
#include "util.fh"
      integer rtdb
c
c     The amplitudes are stored on disk in 2 orderings
c
c     1. Saving amplitudes for restart, C-terms and triples.
c
c     2. B-term (blocked, antisymmtrized SO basis as tij,uv).
c
c     The DRA are 1-dimensional with t1a, t1b, t2aa, t2bb, t2ab, t2ba 
c     stored consecutively.
c
c     Order 1.
c
c     The t1 are stored symmetry packed (as described for an OV pair
c     below).
c
c     The t2 are stored symmetry packed as t(i,a,symi,j,b,symb,symjb).
c     For a given sym(jb)=sym(ia) all amplitudes are contiguous,
c     as are all amplitudes if the symmetry of j is also specified.
c
c     ind = 0
c     do symmetry of jb (= symmetry of ia)
c     .  t(ia,jb,irrep) is a dense matrix
c
c     All OV pairs of a given symmetry are packed as follows
c     (note that orbitals are stored in symmetry order, so that a
c     double loop over irreps and orbitals in the irrep goes over
c     the orbitals in order)
c     
c     Given sym(ia)
c     .  ind = 0
c     .  do symmetry of a
c     .     know symmetry of i and sub-matrix is dense & rectangular
c     .     (#i of symi, #a of syma)
c     .     do a in irrep
c     .        do i in irrep 
c     .           ind = ind + 1
c     .           packed(ind) = unpacked(i,a)
c
c     No index symmetry is exploited so the pure spin amplitudes
c     are stored 4x times too many.  Thus, we will use 3*(O*V)**2 
c     words rather than the mimimum (3/2)*(O*V)**2.
c
c     Reading/writing/accumulating of things is supported via routines
c     which are all collective and transfer to/from a global array.
c     
c     .   t1 amplitudes of given spin into a global array.
c
c     .   t2 amplitudes for given symjb and range of b into a 2-d
c     .   global array of appropriate dimension (below). 
c     
c     To address the files and the corresponding global arrays we have
c     
c     - offsets to the beginning of things in the DRA
c
c     .   file_t1_off(spin)        
c     .   file_t2_off(spini,spina,spinj,spinb,symjb) = offset within file to T2
c     
c     - symmetry packed occupied-virtual pairs are addressed with
c
c     .   ov_off(a,symia,spini,spina) ... 
c     .   and note that ov_off() is dimensioned (max_nmo+1,0:7,2,2)
c     .   so it is OK to compute differences using ov_off(a+1,...)
c
c     .   double precision packed(1:*), unpacked(nocc,nvirt)
c     .   packed(1 + i-o_sym(1,symi,spini) + ov_off(a,symia,spini,spina))
c     .        = unpacked(i,a)
c
c     - dimension information
c     .   ov_len(symia,spini,spina) = length of symmetry packed (i,a) of 
c     .                               given sym and spins
c     .   ov_len_cum(symia,spini,spina) 
c     .       = sum(ov_len(sym,spini,spina),sym=0..symia-1)
c     .                             = ov_off(nmo+1,spini,spina)
c     .   t1_len(spin)              = length of t1 amplitudes of given spin
c     .   t2_len(spini,spinj)       = length of t2 amplitudes of given spin
c     .                               when stored without using index exchange symmetry
c     .   t2_len_min(spini,spinj)   = length of t2 amplitudes of given spin
c     .                               using full index exchange symmetry
c     .   nvariables                = no. of variables in the CCSD
c     .                             = t1_len(1) + t1_len(2) + 
c     .                               t2_len_min(1,1) + t2_len_min(1,2) + 
c     .                               t2_len_min(2,2)
c
c     !!! Note that t2_len* and nvariables might be too large to represent
c     !!! in a 32-bit integer and are declared as LONGINTEGER
c
c     API ... all routines are collective and return true/false on 
c     success/failure.
c
c     uccsdt_init must have been called to initialize the addressing
c     arrays.
c     
c     logical uccsdt_ampfile_init(rtdb)
c
c     logical uccsdt_ampfile_open(filename, d_amp, oscratch)
c     .  Open the file returning the handle. 
c
c     logical uccsdt_ampfile_close(d_amp)
c     .  Close the file.
c
c     logical uccsdt_ampfile_delete(filename)
c     .  Delete the file.
c     
c     logical uccsdt_ampfile_read_t1(d_amp,spin,g_t1,ocreate)
c     .  Reads the t1 amplitudes of specified spin into a global array. 
c     .  If (ocreate) then 
c     .     The global array is created with the minimum dimensions
c     .     which are (nocc,nvirt)
c     .  else
c     .     The global array must exist and meet or exceed the minimum
c     .     dimensions
c     .  endif
c
c     logical uccsdt_ampfile_write_t1(d_amp,spin,g_t1)
c
c     logical uccsdt_ampfile_read_t2(d_amp,spinj,spinb,symjb,blo,bhi,g_t2,ocreate)
c
c     logical uccsdt_ampfile_write_t2(d_amp,spinj,spinb,symjb,blo,bhi,g_t2)
c     
      integer spini, spinj, spina, spinb, spin
      integer symi, syma
      integer symia, symjb
      integer a, ind
      logical oprint
c
      if (dra_init(8, -1, -1, 2*1024*1024) .ne. 0) 
     $     call errquit('uccsdt: initializing disk arrays?',0, DISK_ERR)
      oprint = util_print('fileaddr',print_debug)
c
c     Setup the symmetry packed pair addressing arrays for both pure
c     and mixed spin
c
      call ifill((nw_max_nbf_small+1)*8*2*2,0,ov_len,1)
      do spini = 1, 2
         do spina = 1, 2
            do symia = 0,7
               ind = 0
               do syma = 0, 7
                  symi = ieor(symia,syma)
                  do a = v_sym(1,syma,spina),v_sym(2,syma,spina)
                     ov_off(a,symia,spini,spina) = ind
                     ind = ind + no_sym(symi,spini)
                  end do
               end do
               ov_len(symia,spini,spina) = ind
               ov_off(nmo+1,symia,spini,spina) = ind ! For easy size comp.
            end do
         end do
      end do
      do spini = 1, 2
         do spina = 1, 2
            ov_len_cum(0,spini,spina) = 0
            do symia = 1,7
               ov_len_cum(symia,spini,spina) =  
     $              ov_len_cum(symia-1,spini,spina) + 
     $              ov_len(symia-1,spini,spina)
            end do
         end do
      end do
c      
      if (oprint) then
         write(6,*) ' OV lengths '
         do spini = 1, 2
            do spina = 1, 2
               write(6,71) spini,spina,
     $              (ov_len(symia,spini,spina),symia=0,7)
 71            format(5x,2i5,5x,8i5)
            end do
         end do
         write(6,*) ' OV lengths cummulative '
         do spini = 1, 2
            do spina = 1, 2
               write(6,71) spini,spina,
     $              (ov_len_cum(symia,spini,spina),symia=0,7)
            end do
         end do
         write(6,*) ' OV offsets '
         do spini = 1, 2
            do spina = 1, 2
               do symia = 0,7
                  write(6,72) spini, spina, symia
 72               format(5x,3i5)
                  call ioutput(ov_off(1,symia,spini,spina),
     $                 nc(spina)+no(spina)+1,nmo,1,1,nmo,nmo,1)
               end do
            end do
         end do
      end if
c
c     Size and file offsets of symmetry packed T1 amplitudes
c
      do spin = 1, 2
         t1_len(spin) = 0
         do symia = 0, 7
            t1_len(spin) = t1_len(spin) + ov_len(symia,spin,spin)
         end do
      end do
      file_t1_off(1) = 0
      file_t1_off(2) = t1_len(1)
c
c     Size and file offsets for the T2s
c
c     t2(spini,spina,spinj,spinb,symjb) is a dense 2D array of dimension
c     (ov_len(symia,spini,spina),ov_len(symjb,spinj,spinb)) and begins
c     at file_t2_off(spini,spina,spinj,spinb,symjb)
c
      ind = file_t1_off(2) + t1_len(2)
      do spinb = 1, 2
         do spinj = 1, 2
            do spina = 1, 2
               do spini = 1, 2
                  if ((spini+spinj).eq.(spina+spinb)) then
                     do symjb = 0,7
                        file_t2_off(spini,spina,spinj,spinb,symjb) = ind
                        ind = ind + ov_len(symjb,spini,spina) *
     $                       ov_len(symjb,spinj,spinb)
                     end do
                  else
                     do symjb = 0,7
                        file_t2_off(spini,spina,spinj,spinb,symjb) = 
     $                       -99999999
                     end do
                  end if
               end do
            end do
         end do
      end do
c
c     Lengths of the T2s
c
      do spini = 1, 2
         do spinj = 1, 2
            t2_len(spini,spinj) = 0
            do symjb = 0,7
               t2_len(spini,spinj) = t2_len(spini,spinj) + 
     $              ov_len(symjb,spini,spini) *
     $              ov_len(symjb,spinj,spinj)
            end do
         end do
      end do
c
      if (oprint) then
         write(6,*) ' T1 len ', t1_len
         write(6,*) ' T2 len ', t2_len
         write(6,*) ' T1 file offsets ', file_t1_off
         write(6,*) ' T2 file offsets '
         do spinb = 1, 2
            do spinj = 1, 2
               do spina = 1, 2
                  do spini = 1, 2
                     if ((spini+spinj).eq.(spina+spinb)) then
                        write(6,*) '   ',spini,spina,spinj,spinb
                        do symjb = 0,7
                           write(6,*) '        ',symjb,
     $                          file_t2_off(spini,spina,spinj,
     $                          spinb,symjb)
                        end do
                     end if
                  end do
               end do
            end do
         end do
      end if
c
c      call uccsdt_ampfile_test()
c
      uccsdt_ampfile_init = .true.
c
      end
      logical function uccsdt_ampfile_tidy()
      implicit none
#include "dra.fh"
c
      uccsdt_ampfile_tidy = dra_terminate() .eq. 0
c
      end
      logical function uccsdt_ampfile_close(d_amp, okeep)
      implicit none
#include "dra.fh"
      integer d_amp             ! [input]
      logical okeep             ! [input] True/False for Keep/Delete
c
      if (okeep) then
         uccsdt_ampfile_close = dra_close(d_amp) .eq. 0
      else
         uccsdt_ampfile_close = dra_delete(d_amp) .eq. 0
      endif
c
      end
      logical function uccsdt_ampfile_open(name, filename, mode, d_amp)
      implicit none
#include "cuccsdtP.fh"
#include "dra.fh"
#include "mafdecls.fh"
      character*(*) name        ! [input]
      character*(*) filename    ! [input]
      character*(*) mode        ! [input] new, old
      integer d_amp             ! [output] DRA handle
c
      LONGINTEGER size
c
      if (mode .eq. 'new') then
c
c     Size will be enough for t1a + t1b + t2aaaa + t2bbbb + t2aabb + t2bbaa
c     ... 8/3 times the minimum.
c
         size = t1_len(1) + t1_len(2) + t2_len(1,1) + t2_len(2,2)
     $        +  2*t2_len(1,2) + 2*t2_len(2,1)
c     
         uccsdt_ampfile_open =  dra_create(mt_dbl, size, 1, 
     $        name, filename, DRA_RW, -1, 1, d_amp) .eq. 0
         write(6,*) ' OPENED with size ', size
      else
         uccsdt_ampfile_open = dra_open(filename, DRA_RW, d_amp) .eq. 0
      end if
c
      end
      logical function uccsdt_ampfile_read_t1(d_amp, spin, g_t1, 
     $     ocreate)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "dra.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer d_amp             ! [input]
      integer spin              ! [input]
      integer g_t1              ! [input/output]
      logical ocreate           ! [input]
c
      integer request
c
      if (ocreate) then
         if (.not. ga_create(mt_dbl, t1_len(spin), 1, 'readt1amp',
     $        0, 0, g_t1)) call errquit
     $        ('uccsdt_ampfile_read_t1: ga_create failed?',0, GA_ERR)
      end if
c
      uccsdt_ampfile_read_t1 = dra_read_section(.false., 
     $     g_t1, 1, t1_len(spin), 1, 1,
     $     d_amp, file_t1_off(spin)+1, file_t1_off(spin)+t1_len(spin),
     $     1, 1,
     $     request) .eq. 0
      if (.not. uccsdt_ampfile_read_t1) return
c
      uccsdt_ampfile_read_t1 = dra_wait(request) .eq. 0
c
      end
      logical function uccsdt_ampfile_read_t2(d_amp, 
     $     spini, spina, spinj, spinb, symjb, blo, bhi, 
     $     g_t2, ocreate, distribution)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "dra.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer d_amp             ! [input]
      integer spini, spina, spinj, spinb ! [input]
      integer symjb             ! [input]
      integer blo, bhi          ! [input]
      integer g_t2              ! [input/output]
      logical ocreate           ! [input]
      character*(*) distribution ![input]
c
c     Read the amplitudes for the given spins, symmetries, and 
c     range of b into the global array from the given disk file (d_amp).
c       
c     If (ocreate) then 
c     .   The global array is created with the minimum dimensions
c     .   which are (lenia,lenjb) with
c     .      lenia = ov_len(symjb) ... all OV pairs of given symmetry
c     .      lenjb = ov_off(bhi+1,symjb,spinj,spinb) - 
c     .                ov_off(blo,  symjb,spinj,spinb) ... pairs within range
c     .      if (distribution == 'block') then
c     .         uses the default GA uniform block distribution
c     .      else if (distribution == 'column') then
c     .         distributes only the second index (columns local to a process)
c     .      else
c     .         error
c     .      endif
c     else
c     .   The global array must exist and meet or exceed the minimum 
c     .   dimensions
c     endif
c
      integer request, lenia, lenjb, symia, col, row
      LONGINTEGER offset
c
      symia = symjb
      lenia = ov_len(symia,spini,spina)
      lenjb = ov_off(bhi+1,symjb,spinj,spinb) - 
     $     ov_off(blo,symjb,spinj,spinb)
c
      offset = file_t2_off(spini,spina,spinj,spinb,symjb) +
     $     ov_len(symia,spini,spina) * ov_off(blo,symjb,spinj,spinb)
c
      if (ocreate) then
         if (distribution .eq. 'block') then
            col = -1
            row = -1
         else if (distribution .eq. 'column') then
            col = lenia
            row = -1
         else
            call errquit('uccsdt_ampfile_read_t2: distribution?',0,
     &       INPUT_ERR)
         end if
         if (.not. ga_create(mt_dbl, lenia, lenjb, 'readt2amp',
     $        col, row, g_t2)) call errquit
     $        ('uccsdt_ampfile_read_t2: ga_create failed?',0, GA_ERR)
      end if
c
**      write(6,*) ' OFFSET ', offset+1, offset + lenia*lenjb
      uccsdt_ampfile_read_t2 = dra_read_section(.false., 
     $     g_t2, 1, lenia, 1, lenjb,
     $     d_amp, offset+1, offset + lenia*lenjb, 1,1,
     $     request) .eq. 0
      if (.not. uccsdt_ampfile_read_t2) return
c
      uccsdt_ampfile_read_t2 = dra_wait(request) .eq. 0
c
c      call ga_print(g_t2)
c
      end
      logical function uccsdt_ampfile_write_t1(d_amp, spin, g_t1)
      implicit none
#include "cuccsdtP.fh"
#include "dra.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer d_amp             ! [input]
      integer spin              ! [input]
      integer g_t1              ! [input]
c
      integer request
c
      uccsdt_ampfile_write_t1 = dra_write_section(.false., 
     $     g_t1, 1, t1_len(spin), 1, 1,
     $     d_amp, file_t1_off(spin)+1, file_t1_off(spin)+t1_len(spin),
     $     1, 1,
     $     request) .eq. 0
      if (.not. uccsdt_ampfile_write_t1) return
c
      uccsdt_ampfile_write_t1 = dra_wait(request) .eq. 0
c
      end
      logical function uccsdt_ampfile_write_t2(d_amp, 
     $     spini, spina, spinj, spinb, symjb, blo, bhi, 
     $     g_t2)
      implicit none
#include "cuccsdtP.fh"
#include "dra.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer d_amp             ! [input]
      integer spini, spina, spinj, spinb ! [input]
      integer symjb             ! [input]
      integer blo, bhi          ! [input]
      integer g_t2              ! [input]
c
c     Write the amplitudes for the given spins, symmetries, and 
c     range of b into the global array from the given disk file (d_amp).
c       
      integer request, lenia, lenjb, symia
      LONGINTEGER offset
c
      symia = symjb
      lenia = ov_len(symia,spini,spina)
      lenjb = ov_off(bhi+1,symjb,spinj,spinb) - 
     $     ov_off(blo,symjb,spinj,spinb)
c
      offset = file_t2_off(spini,spina,spinj,spinb,symjb) +
     $     ov_len(symia,spini,spina) * ov_off(blo,symjb,spinj,spinb)
c
      uccsdt_ampfile_write_t2 = dra_write_section(.false., 
     $     g_t2, 1, lenia, 1, lenjb,
     $     d_amp, 
     $     offset + 1, offset + lenia*lenjb, 1,1,
     $     request) .eq. 0
      if (.not. uccsdt_ampfile_write_t2) return
c
      uccsdt_ampfile_write_t2 = dra_wait(request) .eq. 0
c
      end
      subroutine uccsdt_ampfile_test()
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
c
      integer d_amp, g_t1, g_t2
      integer symia, symjb
      integer symi, syma, symj, symb, spini, spina, spinj, spinb
      integer blo, bhi, ia, jb, i, a, j, b
      double precision test, test2
c
      if (.not. uccsdt_ampfile_open('fred','/tmp/fred','new',d_amp))
     $     call errquit('amplitudes: open failed',0, DISK_ERR)
c
      symjb = 1
      do spinb = 1, 2
         spinj = spinb
         if (.not. uccsdt_ampfile_read_t1(d_amp, spinb, g_t1, 
     $        .true.)) call errquit('read of t1 failed',0, DISK_ERR)
         jb = 0
         do symb = 0,7
            symj = ieor(symb,symjb)
            do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
               do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
                  jb = jb + 1
                  test = j + 100*b
                  call ga_put(g_t1, jb, jb, 1, 1, test, 1)
               end do
            end do
         end do
         if (.not. uccsdt_ampfile_write_t1(d_amp, spinb, g_t1)) 
     $        call errquit('write of t1 failed',0, DISK_ERR)
         if (.not. ga_destroy(g_t1)) call errquit('GA?',0, GA_ERR)
      end do
c
      do spinb = 1, 2
         do spinj = 1, 2
            do spina = 1, 2
               do spini = 1, 2
                  if ((spini+spinj).eq.(spina+spinb)) then
                     write(6,*) ' SPINS ', spini,spina,spinj,spinb
                     do symjb = 0, 7
                        symia = symjb
                        write(6,*) '        symjb', symjb, 
     $                       ov_len(symjb,spinj,spinb),
     $                       ov_len(symia,spini,spina)
                        if ( ov_len(symjb,spinj,spinb).gt.0 .and.
     $                       ov_len(symia,spini,spina).gt.0) then
                        blo = nc(spinb)+no(spinb)+1
                        bhi = nmo
                        if (.not. uccsdt_ampfile_read_t2(d_amp,
     $                       spini, spina, spinj, spinb, symjb,
     $                       blo, bhi, g_t2, .true., 'block')) 
     $                       call errquit
     $                       ('read of all failed',0, DISK_ERR)
                        jb = 0
                        do symb = 0,7
                           symj = ieor(symb,symjb)
                           do b = v_sym(1,symb,spinb),
     $                          v_sym(2,symb,spinb)
                              if (jb.ne.ov_off(b,symjb,spinj,spinb))
     $                             call errquit('bad jb',0, INPUT_ERR)
                              do j = o_sym(1,symj,spinj),
     $                             o_sym(2,symj,spinj)
                                 jb = jb + 1
c
                                 ia = 0
                                 do syma = 0,7
                                    symi = ieor(syma,symia)
                                    do a = v_sym(1,syma,spina),
     $                                   v_sym(2,syma,spina)
                                       if (ia.ne.ov_off(a,symia,spini,
     $                                      spina))
     $                                      call errquit('bad ia',0,
     &       INPUT_ERR)
                                       do i = o_sym(1,symi,spini),
     $                                      o_sym(2,symi,spini)
                                          ia = ia + 1
c
                                          test = i+100*(a+100*(j+100*b))
                                          call ga_put(g_t2,
     $                                         ia,ia,jb,jb,test,1)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                        if (.not. uccsdt_ampfile_write_t2(d_amp,
     $                       spini, spina, spinj, spinb, symjb,
     $                       blo, bhi, g_t2)) call errquit
     $                       ('write of all failed',0, DISK_ERR)
                        if (.not. ga_destroy(g_t2)) call errquit
     $                       ('ga?',0, GA_ERR)
                        endif
                     end do
                  end if
               end do
            end do
         end do
      end do
c
c     All of the amplitudes have now been written out ... read
c     them back in a different order and in smaller pieces.
c
      symjb = 1
      do spinb = 1, 2
         spinj = spinb
         if (.not. uccsdt_ampfile_read_t1(d_amp, spinb, g_t1, 
     $        .true.)) call errquit('read of t1 failed',0, DISK_ERR)
         jb = 0
         do symb = 0,7
            symj = ieor(symb,symjb)
            do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
               do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
                  jb = jb + 1
                  test = j + 100*b
                  call ga_get(g_t1, jb, jb, 1, 1, test2, 1)
                  if (test2.ne.test) call errquit('wrong t1',0,
     &       INPUT_ERR)
               end do
            end do
         end do
         if (.not. ga_destroy(g_t1)) call errquit('GA?',0, GA_ERR)
      end do
c
      do spini = 1, 2
         do spinb = 1, 2
            do spina = 1, 2
               do spinj = 1, 2
                  if ((spini+spinj).eq.(spina+spinb)) then
                     write(6,*) ' SPINS ', spini,spina,spinj,spinb
                     do symjb = 7,0,-1
                        symia = symjb
                        write(6,*) '        symjb', symjb, 
     $                       ov_len(symjb,spinj,spinb),
     $                       ov_len(symia,spini,spina)
                        if ( ov_len(symjb,spinj,spinb).gt.0 .and.
     $                       ov_len(symia,spini,spina).gt.0) then
                        do symb = 0,7
                           symj = ieor(symb,symjb)
                           if ( nv_sym(symb,spinb)*
     $                          no_sym(symj,spinj) .gt. 0) then
                           blo = v_sym(1,symb,spinb)
                           bhi = v_sym(2,symb,spinb)
                           if (.not. uccsdt_ampfile_read_t2(d_amp,
     $                          spini, spina, spinj, spinb, symjb,
     $                          blo, bhi, g_t2, .true.,'block')) 
     $                          call errquit
     $                          ('read of all failed',0, DISK_ERR)
                           jb = 0
                           do b = v_sym(1,symb,spinb),
     $                          v_sym(2,symb,spinb)
                              if (jb.ne.ov_off(b,symjb,spinj,spinb)-
     $                             ov_off(blo,symjb,spinj,spinb))
     $                             call errquit('bad jb2 ',0, INPUT_ERR)
                              do j = o_sym(1,symj,spinj),
     $                             o_sym(2,symj,spinj)
                                 jb = jb + 1
c
                                 ia = 0
                                 do syma = 0,7
                                    symi = ieor(syma,symia)
                                    do a = v_sym(1,syma,spina),
     $                                   v_sym(2,syma,spina)
                                       if (ia.ne.ov_off(a,symia,spini,
     $                                      spina))
     $                                      call errquit('bad ia',0,
     &       INPUT_ERR)
                                       do i = o_sym(1,symi,spini),
     $                                      o_sym(2,symi,spini)
                                          ia = ia + 1
c
                                          test = i+100*(a+100*(j+100*b))
                                          call ga_get(g_t2,
     $                                         ia,ia,jb,jb,test2,1)
                                          if (test.ne.test2) call 
     $                                         errquit('wrong!',0,
     &       INPUT_ERR)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                           if (.not. ga_destroy(g_t2)) call errquit
     $                          ('ga2?',0, GA_ERR)
                           end if
                        end do
                        endif
                     end do
                  end if
               end do
            end do
         end do
      end do
c
      if (.not. uccsdt_ampfile_close(d_amp,.false.)) call errquit
     $     ('close delete failed',0, DISK_ERR)
c
      end
      subroutine uccsdt_ampfile_test_write_t1(t1a, t1b,
     &                                  oloa, ohia, vloa, vhia, 
     &                                  olob, ohib, vlob, vhib,
     &                                  d_amp)
      implicit none
#include "errquit.fh"
c
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
      integer d_amp, g_t1
      integer symjb
      integer symj, symb, spinj, spinb
      integer jb, j, b
      integer oloa, ohia, olob, ohib, vloa, vhia, vlob, vhib
      integer max_len_t1buf, l_t1buf, k_t1buf
      double precision t1a(vloa:vhia,oloa:ohia)
      double precision t1b(vlob:vhib,olob:ohib)
      double precision test
c
c     write t1a, t1b
c
c      call debug_print_t1('k_t1a', t1a, vloa, vhia, oloa, ohia)
c      call debug_print_t1('k_t1b', t1b, vlob, vhib, olob, ohib)
c
c
c     putting t1 from ma -> ga -> dra (for testing purposes)
c
c     Get a chunk of MA (max(len_t1(1),len_t1(2)))
c
      max_len_t1buf = max(t1_len(1),t1_len(2))
      if (.not. ma_push_get(mt_dbl, max_len_t1buf, 't buffer',
     $     l_t1buf, k_t1buf)) call errquit(' t1 read buffer: ma?', 1,
     &       MA_ERR)
      symjb = 0
      do spinb = 1, 2
         spinj = spinb
         if (.not. uccsdt_ampfile_read_t1(d_amp, spinb, g_t1, 
     $        .true.)) call errquit('read of t1 failed',0, DISK_ERR)
         jb = 0
c         do symb = 0,7
c            symj = ieor(symb,symjb)
c            write(6,*)' spinb, spinj: ',
c     &                  spinb, spinj
c            write(6,*)' symb, symjb, symj: ',
c     &                  symb, symjb, symj
c            write(6,*)' v_sym(1,symb,spinb), v_sym(2,symb,spinb): ',
c     &                  v_sym(1,symb,spinb), v_sym(2,symb,spinb)
c            write(6,*)' o_sym(1,symj,spinj),o_sym(2,symj,spinj): ',
c     &                  o_sym(1,symj,spinj),o_sym(2,symj,spinj)
c            do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
c               do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
c                  jb = jb + 1
c                  if (spinb.eq.1)then
c                     test = t1a(b,j)
c                  elseif (spinb.eq.2)then
c                     test = t1b(b,j)
c                  endif
c                  write(6,*)' b, j, jb, test: ',b,j,jb,test
c            
c                  call ga_put(g_t1, jb, jb, 1, 1, test, 1)
c               end do
c            end do
c         end do
         if (spinb.eq.1) then
           write(6,*)' packing t1a buffer '
           call pack_t1buf(t1a, dbl_mb(k_t1buf),
     &          t1_len(1), symjb, jb, 
     &          oloa, ohia, vloa, vhia, 
     &          spinb, spinj)
         elseif (spinb.eq.2) then
           write(6,*)' packing t1b buffer '
           call pack_t1buf(t1b, dbl_mb(k_t1buf),
     &          t1_len(2), symjb, jb, 
     &          olob, ohib, vlob, vhib, 
     &          spinb, spinj)
         endif
         call ga_put(g_t1, 1, jb, 1, 1, dbl_mb(k_t1buf), jb)
         if (.not. uccsdt_ampfile_write_t1(d_amp, spinb, g_t1)) 
     $        call errquit('write of t1 failed',0, GA_ERR)
c
c         call ga_print(g_t1)
c
         if (.not. ga_destroy(g_t1)) call errquit('GA?',0, GA_ERR)
      end do
c
      if (.not. ma_pop_stack(l_t1buf)) call errquit('t1buf?',0, MA_ERR)
c
      return
      end
      subroutine uccsdt_ampfile_test_write_t2(t2aa, t2bb, t2ab,
     &                                  oloa, ohia, vloa, vhia, 
     &                                  olob, ohib, vlob, vhib,
     &                                  d_amp)
      implicit none
#include "errquit.fh"
c
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
      integer d_amp, g_t2
      integer symia, symjb
      integer blo, bhi, ia, jb, i, a, j, b
      integer symi, syma, symj, symb, spini, spina, spinj, spinb
      integer oloa, ohia, olob, ohib, vloa, vhia, vlob, vhib
      integer lenia, lenjb, max_lenjb, l_t2buf, k_t2buf, spinsum
      double precision t2aa(vloa:vhia,vloa:vhia,oloa:ohia,oloa:ohia)
      double precision t2bb(vlob:vhib,vlob:vhib,olob:ohib,olob:ohib)
      double precision t2ab(vloa:vhia,vlob:vhib,oloa:ohia,olob:ohib)
      double precision test
c
c     write t2aa, t2bb, t2ab
c
c      call debug_print_t2('k_t2aa', t2aa, vloa, vhia, vloa, vhia, 
c     &                                    oloa, ohia, oloa, ohia)
c      call debug_print_t2('k_t2bb', t2bb, vlob, vhib, vlob, vhib, 
c     &                                    olob, ohib, olob, ohib)
c      call debug_print_t2('k_t2ab', t2ab, vloa, vhia, vlob, vhib, 
c     &                                    oloa, ohia, olob, ohib)
c
c     putting t2 from ma -> ga -> dra (for testing purposes)
c     Get a chunk of MA (needs to be max_lenia*max_lenjb (for now))
c
      max_lenjb = 0
      do spinb = 1, 2
         do spinj = 1, 2
           do symjb = 0, 7
             max_lenjb = max(ov_len(symjb,spinj,spinb),max_lenjb)
           enddo
         enddo
      enddo
      if (.not. ma_push_get(mt_dbl, max_lenjb*max_lenjb, 't buffer',
     $     l_t2buf, k_t2buf)) call errquit(' t2 read buffer: ma?', 1,
     &       MA_ERR)
c      write(6,*)' max_lenjb = ', max_lenjb
c
      do spinb = 1, 2
         do spinj = 1, 2
            do spina = 1, 2
               do spini = 1, 2
                  if ((spini+spinj).eq.(spina+spinb)) then
                     spinsum = spinb + spinj + spina + spini
                     write(6,*) ' SPINS ', spini,spina,spinj,spinb
                     do symjb = 0, 7
                        symia = symjb
c                        write(6,*) '        symjb', symjb, 
c     $                       ov_len(symjb,spinj,spinb),
c     $                       ov_len(symia,spini,spina)
                        if ( ov_len(symjb,spinj,spinb).gt.0 .and.
     $                       ov_len(symia,spini,spina).gt.0) then
                        blo = nc(spinb)+no(spinb)+1
                        bhi = nmo
                        if (.not. uccsdt_ampfile_read_t2(d_amp,
     $                       spini, spina, spinj, spinb, symjb,
     $                       blo, bhi, g_t2, .true., 'block')) 
     $                       call errquit
     $                       ('read of all failed',0, DISK_ERR)
c                        jb = 0
c                        do symb = 0,7
c                           symj = ieor(symb,symjb)
c                           do b = v_sym(1,symb,spinb),
c     $                          v_sym(2,symb,spinb)
c                              if (jb.ne.ov_off(b,symjb,spinj,spinb))
c     $                             call errquit('bad jb',0)
c                              do j = o_sym(1,symj,spinj),
c     $                             o_sym(2,symj,spinj)
c                                 jb = jb + 1
cc
c                                 ia = 0
c                                 do syma = 0,7
c                                    symi = ieor(syma,symia)
c                                    do a = v_sym(1,syma,spina),
c     $                                   v_sym(2,syma,spina)
c                                       if (ia.ne.ov_off(a,symia,spini,
c     $                                      spina))
c     $                                      call errquit('bad ia',0)
c                                       do i = o_sym(1,symi,spini),
c     $                                      o_sym(2,symi,spini)
c                                          ia = ia + 1
c                                          if (spinsum.eq.4)then
c                                             test = t2aa(a,b,i,j)
c                                          elseif (spinsum.eq.8)then
c                                             test = t2bb(a,b,i,j)
c                                          else
c                                             test = t2ab(a,b,i,j)
c                                             if (spina.ne.spini)
c     &                                          test = -test
c                                          endif
cc
c                write(6,*)' spinsum, ia, jb, a, b, i, j, test ',
c     &                      spinsum, ia, jb, a, b, i, j, test
c                                          call ga_put(g_t2,
c     $                                         ia,ia,jb,jb,test,1)
c                                       end do
c                                    end do
c                                 end do
c                              end do
c                           end do
c                        end do
c
c                       stuff t2 elements into a buffer to be put into ga
c
                        lenia = ov_len(symia,spini,spina)
                        lenjb = ov_off(bhi+1,symjb,spinj,spinb) - 
     $                          ov_off(blo,symjb,spinj,spinb)
c                        write(6,*)' putting lenia, lenjb ',
c     &                                      lenia, lenjb
                        if (spinsum.eq.4)then
                          write(6,*)' packing t2aa buffer '
                          call pack_t2buf(t2aa, dbl_mb(k_t2buf),
     &                         lenia, lenjb, symjb,
     &                         oloa, ohia, vloa, vhia, 
     &                         oloa, ohia, vloa, vhia,
     &                         spina, spinb, spini, spinj,.false.)
                        elseif (spinsum.eq.8)then
                          write(6,*)' packing t2bb buffer '
                          call pack_t2buf(t2bb, dbl_mb(k_t2buf),
     &                         lenia, lenjb, symjb,
     &                         olob, ohib, vlob, vhib, 
     &                         olob, ohib, vlob, vhib,
     &                         spina, spinb, spini, spinj,.false.)
                        else
                          write(6,*)' packing t2ab buffer '
                          if ((spinj.gt.spini).and.
     &                        (spinb.gt.spina))then
                            call pack_t2buf(t2ab, dbl_mb(k_t2buf),
     &                           lenia, lenjb, symjb,
     &                           oloa, ohia, vloa, vhia, 
     &                           olob, ohib, vlob, vhib,
     &                           spina, spinb, spini, spinj,.false.)
                          elseif ((spinj.lt.spini).and.
     &                           (spinb.lt.spina))then
                            call pack_t2buf(t2ab, dbl_mb(k_t2buf),
     &                           lenia, lenjb, symjb,
     &                           oloa, ohia, vloa, vhia, 
     &                           olob, ohib, vlob, vhib,
     &                           spina, spinb, spini, spinj,.true.)
                          endif
                        endif
c
                        call ga_put(g_t2, 1, lenia, 1, lenjb, 
     &                              dbl_mb(k_t2buf), lenia)
c
                        if (.not. uccsdt_ampfile_write_t2(d_amp,
     $                       spini, spina, spinj, spinb, symjb,
     $                       blo, bhi, g_t2)) call errquit
     $                       ('write of all failed',0, DISK_ERR)
c
                        if (.not. ga_destroy(g_t2)) call errquit
     $                       ('ga?',0, GA_ERR)
                        endif
                     end do
                  end if
               end do
            end do
         end do
      end do
c
      if (.not. ma_pop_stack(l_t2buf)) call errquit('t2buf?',0, MA_ERR)
c
      return
      end
      subroutine debug_print_t2(string, t2, vlo1, vhi1, vlo2, vhi2, 
     &                          olo1, ohi1, olo2, ohi2)
      implicit none
      character*(*) string
      integer vlo1, vhi1, vlo2, vhi2, olo1, ohi1, olo2, ohi2
      double precision t2(vlo1:vhi1,vlo2:vhi2,olo1:ohi1,olo2:ohi2)
c
      integer a, b, i, j
      write(6,*) ' DEBUG FOR ', string, vlo1, vhi1, vlo2, vhi2, 
     &                          olo1, ohi1, olo2, ohi2
      do j = olo2, ohi2
         do i = olo1, ohi1
            do b = vlo2, vhi2
               do a = vlo1, vhi1
                  if (abs(t2(a,b,i,j)).gt.1e-6) then
                     write(6,1) a,b,i,j,t2(a,b,i,j)
 1                   format(4i5,2x,f14.8)
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine debug_print_t1(string, t1, vlo1, vhi1, olo1, ohi1)
      implicit none
      character*(*) string
      integer vlo1, vhi1, olo1, ohi1
      double precision t1(vlo1:vhi1,olo1:ohi1)
c
      integer a, i
      write(6,*) ' DEBUG FOR ', string, vlo1, vhi1, olo1, ohi1
      do i = olo1, ohi1
        do a = vlo1, vhi1
          if (abs(t1(a,i)).gt.1e-6) then
            write(6,1) a,i,t1(a,i)
    1       format(2i5,2x,f14.8)
          end if
        end do
      end do
c
      end
      subroutine amps_copy(t1a, t1b,
     &               t2aa, t2bb, t2ab,
     &               t1ar, t1br,
     &               t2aar, t2bbr, t2abr,
     &               nva, nvb, noa, nob, nbf,
     &               oloa, ohia, vloa, vhia, 
     &               olob, ohib, vlob, vhib)
      implicit none
      integer nva, nvb, noa, nob, nbf,
     &        oloa, ohia, vloa, vhia, 
     &        olob, ohib, vlob, vhib
      double precision t1a(nva,noa)
      double precision t1b(nvb,nob)
      double precision t2aa(nbf,nbf,nbf,nbf)
      double precision t2bb(nbf,nbf,nbf,nbf)
      double precision t2ab(nbf,nbf,nbf,nbf)
      double precision t1ar(vloa:vhia,oloa:ohia)
      double precision t1br(vlob:vhib,olob:ohib)
      double precision t2aar(vloa:vhia,vloa:vhia,oloa:ohia,oloa:ohia)
      double precision t2bbr(vlob:vhib,vlob:vhib,olob:ohib,olob:ohib)
      double precision t2abr(vloa:vhia,vlob:vhib,oloa:ohia,olob:ohib)
      integer v1, v2, o1, o2
c
      do o2 = oloa, ohia
        do o1 = oloa, ohia
          do v2 = vloa, vhia
            do v1 = vloa, vhia
              t2aar(v1,v2,o1,o2) = t2aa(v1,v2,o1,o2)
            enddo
          enddo
        enddo
      enddo
      do o2 = olob, ohib
        do o1 = olob, ohib
          do v2 = vlob, vhib
            do v1 = vlob, vhib
              t2bbr(v1,v2,o1,o2) = t2bb(v1,v2,o1,o2)
            enddo
          enddo
        enddo
      enddo
      do o2 = olob, ohib
        do o1 = oloa, ohia
          do v2 = vlob, vhib
            do v1 = vloa, vhia
              t2abr(v1,v2,o1,o2) = t2ab(v1,v2,o1,o2)
            enddo
          enddo
        enddo
      enddo
      do o1 = oloa, ohia
        do v1 = vloa, vhia
          t1ar(v1,o1) = t1a(v1-noa,o1)
        enddo
      enddo
      do o1 = olob, ohib
        do v1 = vlob, vhib
          t1br(v1,o1) = t1b(v1-nob,o1)
        enddo
      enddo
      return
      end
      subroutine amps_diff(t1a, t1b,
     &               t2aa, t2bb, t2ab,
     &               t1ar, t1br,
     &               t2aar, t2bbr, t2abr,
     &               nva, nvb, noa, nob, nbf,
     &               oloa, ohia, vloa, vhia, 
     &               olob, ohib, vlob, vhib)
      implicit none
      integer nva, nvb, noa, nob, nbf,
     &        oloa, ohia, vloa, vhia, 
     &        olob, ohib, vlob, vhib
      double precision t1a(nva,noa)
      double precision t1b(nvb,nob)
      double precision t2aa(nbf,nbf,nbf,nbf)
      double precision t2bb(nbf,nbf,nbf,nbf)
      double precision t2ab(nbf,nbf,nbf,nbf)
      double precision t1ar(vloa:vhia,oloa:ohia)
      double precision t1br(vlob:vhib,olob:ohib)
      double precision t2aar(vloa:vhia,vloa:vhia,oloa:ohia,oloa:ohia)
      double precision t2bbr(vlob:vhib,vlob:vhib,olob:ohib,olob:ohib)
      double precision t2abr(vloa:vhia,vlob:vhib,oloa:ohia,olob:ohib)
      integer v1, v2, o1, o2
c
c     diff 2 sets of t1 and t2 amplitudes to make sure they are the same
c
      write(6,*)' diffing 2 sets of amplitudes '
      do o2 = oloa, ohia
        do o1 = oloa, ohia
          do v2 = vloa, vhia
            do v1 = vloa, vhia
              if (dabs(t2aar(v1,v2,o1,o2)-t2aa(v1,v2,o1,o2)).gt.1d-8)
     &             write(6,*)' v1, v2, o1, o2, t2aar, t2aa: ',
     &          v1, v2, o1, o2, t2aar(v1,v2,o1,o2), t2aa(v1,v2,o1,o2)
            enddo
          enddo
        enddo
      enddo
      do o2 = olob, ohib
        do o1 = olob, ohib
          do v2 = vlob, vhib
            do v1 = vlob, vhib
              if (dabs(t2bbr(v1,v2,o1,o2)-t2bb(v1,v2,o1,o2)).gt.1d-8)
     &             write(6,*)' v1, v2, o1, o2, t2bbr, t2bb: ',
     &          v1, v2, o1, o2, t2bbr(v1,v2,o1,o2), t2bb(v1,v2,o1,o2)
            enddo
          enddo
        enddo
      enddo
      do o2 = olob, ohib
        do o1 = oloa, ohia
          do v2 = vlob, vhib
            do v1 = vloa, vhia
              if (dabs(t2abr(v1,v2,o1,o2)-t2ab(v1,v2,o1,o2)).gt.1d-8)
     &             write(6,*)' v1, v2, o1, o2, t2abr, t2ab: ',
     &          v1, v2, o1, o2, t2abr(v1,v2,o1,o2), t2ab(v1,v2,o1,o2)
            enddo
          enddo
        enddo
      enddo
      do o1 = oloa, ohia
        do v1 = vloa, vhia
          t1ar(v1,o1) = t1a(v1-noa,o1)
          if (dabs(t1ar(v1,o1)-t1a(v1-noa,o1)).gt.1d-8)
     &         write(6,*)' v1, o1, t1ar, t1a: ',
     &         v1, o1, t1ar(v1,o1), t1a(v1-noa,o1)
        enddo
      enddo
      do o1 = olob, ohib
        do v1 = vlob, vhib
          if (dabs(t1br(v1,o1)-t1b(v1-nob,o1)).gt.1d-8)
     &         write(6,*)' v1, o1, t1br, t1b: ',
     &         v1, o1, t1br(v1,o1), t1b(v1-nob,o1)
        enddo
      enddo
      return
      end
      subroutine uccsdt_ampfile_test_read_t1(t1a, t1b,
     &                                  oloa, ohia, vloa, vhia, 
     &                                  olob, ohib, vlob, vhib,
     &                                  d_amp)
      implicit none
c
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
      integer d_amp
      integer oloa, ohia, olob, ohib, vloa, vhia, vlob, vhib
      double precision t1a(vloa:vhia,oloa:ohia)
      double precision t1b(vlob:vhib,olob:ohib)
c
      call uccsdt_get_unpacked_t1(d_amp, 1, t1a)
      call uccsdt_get_unpacked_t1(d_amp, 2, t1b)
c
c     write t1a, t1b
c
c      call debug_print_t1('k_t1a', t1a, vloa, vhia, oloa, ohia)
c      call debug_print_t1('k_t1b', t1b, vlob, vhib, olob, ohib)
      return
      end
      subroutine uccsdt_ampfile_test_read_t2(t2aa, t2bb, t2ab,
     &                                  oloa, ohia, vloa, vhia, 
     &                                  olob, ohib, vlob, vhib,
     &                                  d_amp)
      implicit none
#include "errquit.fh"
c
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
      integer d_amp, g_t2
      integer symia, symjb
      integer blo, bhi, ia, jb, i, a, j, b
      integer symi, syma, symj, symb, spini, spina, spinj, spinb
      integer oloa, ohia, olob, ohib, vloa, vhia, vlob, vhib
      integer lenia, lenjb, max_lenjb, l_t2buf, k_t2buf, spinsum
      double precision t2aa(vloa:vhia,vloa:vhia,oloa:ohia,oloa:ohia)
      double precision t2bb(vlob:vhib,vlob:vhib,olob:ohib,olob:ohib)
      double precision t2ab(vloa:vhia,vlob:vhib,oloa:ohia,olob:ohib)
      double precision test
c
c     putting t2 from dra -> ga -> ma (for testing purposes)
c     Get a chunk of MA (needs to be max_lenia*max_lenjb (for now))
c
      max_lenjb = 0
      do spinb = 1, 2
         do spinj = 1, 2
           do symjb = 0, 7
             max_lenjb = max(ov_len(symjb,spinj,spinb),max_lenjb)
           enddo
         enddo
      enddo
      if (.not. ma_push_get(mt_dbl, max_lenjb*max_lenjb, 't buffer',
     $     l_t2buf, k_t2buf)) call errquit(' t2 read buffer: ma?', 1,
     &       MA_ERR)
      call dfill (max_lenjb*max_lenjb, 0.0d0, dbl_mb(k_t2buf),1)
c      write(6,*)' max_lenjb = ', max_lenjb
      do spinb = 1, 2
         do spinj = 1, 2
            do spina = 1, 2
               do spini = 1, 2
                  if ((spini+spinj).eq.(spina+spinb)) then
                     spinsum = spinb + spinj + spina + spini
c                     write(6,*) ' SPINS ', spini,spina,spinj,spinb
                     do symjb = 0, 7
                        symia = symjb
c                        write(6,*) '        symjb', symjb, 
c     $                       ov_len(symjb,spinj,spinb),
c     $                       ov_len(symia,spini,spina)
                        if ( ov_len(symjb,spinj,spinb).gt.0 .and.
     $                       ov_len(symia,spini,spina).gt.0) then
                        blo = nc(spinb)+no(spinb)+1
                        bhi = nmo
                        if (.not. uccsdt_ampfile_read_t2(d_amp,
     $                       spini, spina, spinj, spinb, symjb,
     $                       blo, bhi, g_t2, .true., 'block')) 
     $                       call errquit
     $                       ('read of t2 block failed',0, DISK_ERR)
c
c                        call ga_print(g_t2)
c
                        lenia = ov_len(symia,spini,spina)
                        lenjb = ov_off(bhi+1,symjb,spinj,spinb) - 
     $                          ov_off(blo,symjb,spinj,spinb)
c                        write(6,*)' getting lenia, lenjb ',
c     &                                      lenia, lenjb
                        call ga_get(g_t2, 1, lenia, 1, lenjb, 
     &                              dbl_mb(k_t2buf), lenia)
c
c                       stuff these elements into t2
c
                        if (spinsum.eq.4)then
                          write(6,*)' unpacking t2aa buffer '
                          call unpack_t2buf(t2aa, dbl_mb(k_t2buf),
     &                         lenia, lenjb, symjb,
     &                         oloa, ohia, vloa, vhia, 
     &                         oloa, ohia, vloa, vhia,
     &                         spina, spinb, spini, spinj,.false.)
                        elseif (spinsum.eq.8)then
                          write(6,*)' unpacking t2bb buffer '
                          call unpack_t2buf(t2bb, dbl_mb(k_t2buf),
     &                         lenia, lenjb, symjb,
     &                         olob, ohib, vlob, vhib, 
     &                         olob, ohib, vlob, vhib,
     &                         spina, spinb, spini, spinj,.false.)
                        else
                          write(6,*)' unpacking t2ab buffer '
                          if ((spinj.gt.spini).and.
     &                        (spinb.gt.spina))then
                            call unpack_t2buf(t2ab, dbl_mb(k_t2buf),
     &                           lenia, lenjb, symjb,
     &                           oloa, ohia, vloa, vhia, 
     &                           olob, ohib, vlob, vhib,
     &                           spina, spinb, spini, spinj,.false.)
                          elseif ((spinj.lt.spini).and.
     &                           (spinb.lt.spina))then
                            call unpack_t2buf(t2ab, dbl_mb(k_t2buf),
     &                           lenia, lenjb, symjb,
     &                           oloa, ohia, vloa, vhia, 
     &                           olob, ohib, vlob, vhib,
     &                           spina, spinb, spini, spinj,.true.)
                          endif
                        endif
c
                        if (.not. ga_destroy(g_t2)) call errquit
     $                       ('ga?',0, GA_ERR)
                        endif
                     end do
                  end if
                  
               end do
            end do
         end do
      end do
c
      if (.not. ma_pop_stack(l_t2buf)) call errquit('t2buf?',0, MA_ERR)
c
c     write t2aa, t2bb, t2ab
c
c      call debug_print_t2('dra_t2aa', t2aa, vloa, vhia, vloa, vhia, 
c     &                                    oloa, ohia, oloa, ohia)
c      call debug_print_t2('dra_t2bb', t2bb, vlob, vhib, vlob, vhib, 
c     &                                    olob, ohib, olob, ohib)
c      call debug_print_t2('dra_t2ab', t2ab, vloa, vhia, vlob, vhib, 
c     &                                    oloa, ohia, olob, ohib)
      end
      subroutine unpack_t2buf(t2, t2buf,
     &                        lenia, lenjb, symjb,
     &                        olo1, ohi1, vlo1, vhi1, 
     &                        olo2, ohi2, vlo2, vhi2,
     &                        spina, spinb, spini, spinj, oflip)
      implicit none
#include "errquit.fh"
c
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
      integer lenia, lenjb
      integer symia, symjb
      integer ia, jb, i, a, j, b
      integer symi, syma, symj, symb, spini, spina, spinj, spinb
      integer olo1, ohi1, olo2, ohi2, vlo1, vhi1, vlo2, vhi2
      double precision t2(vlo1:vhi1,vlo2:vhi2,olo1:ohi1,olo2:ohi2)
      double precision t2buf(lenia,lenjb)
      logical oflip
c      call output(t2buf, 1, lenia, 1, lenjb, lenia, lenjb, 1)
      symia = symjb
      jb = 0
      do symb = 0,7
        symj = ieor(symb,symjb)
        do b = v_sym(1,symb,spinb),
     $       v_sym(2,symb,spinb)
          if (jb.ne.ov_off(b,symjb,spinj,spinb))
     $         call errquit('bad jb',0, INPUT_ERR)
          do j = o_sym(1,symj,spinj),
     $         o_sym(2,symj,spinj)
            jb = jb + 1
c           
            ia = 0
            do syma = 0,7
              symi = ieor(syma,symia)
              do a = v_sym(1,syma,spina),
     $             v_sym(2,syma,spina)
                if (ia.ne.ov_off(a,symia,spini,
     $               spina))
     $               call errquit('bad ia',0, INPUT_ERR)
                do i = o_sym(1,symi,spini),
     $               o_sym(2,symi,spini)
                  ia = ia + 1
c                  write(6,*)' ia, jb, a, b, i, j, t2buf(ia,jb) ',
c     &                        ia, jb, a, b, i, j, t2buf(ia,jb)
                  if(oflip)then
                    t2(b,a,j,i) = t2buf(ia,jb)
                  else
                    t2(a,b,i,j) = t2buf(ia,jb)
                  endif
c
                end do
              end do
            end do
          end do
        end do
      end do
      return
      end
      subroutine pack_t2buf(t2, t2buf,
     &                      lenia, lenjb, symjb,
     &                      olo1, ohi1, vlo1, vhi1, 
     &                      olo2, ohi2, vlo2, vhi2,
     &                      spina, spinb, spini, spinj, oflip)
      implicit none
#include "errquit.fh"
c
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
      integer lenia, lenjb
      integer symia, symjb
      integer ia, jb, i, a, j, b
      integer symi, syma, symj, symb, spini, spina, spinj, spinb
      integer olo1, ohi1, olo2, ohi2, vlo1, vhi1, vlo2, vhi2
      double precision t2(vlo1:vhi1,vlo2:vhi2,olo1:ohi1,olo2:ohi2)
      double precision t2buf(lenia,lenjb)
      logical oflip
      call dfill (lenia*lenjb, 0.0d0, t2buf,1)
      symia = symjb
      jb = 0
      do symb = 0,7
        symj = ieor(symb,symjb)
        do b = v_sym(1,symb,spinb),
     $       v_sym(2,symb,spinb)
          if (jb.ne.ov_off(b,symjb,spinj,spinb))
     $         call errquit('bad jb',0, INPUT_ERR)
          do j = o_sym(1,symj,spinj),
     $         o_sym(2,symj,spinj)
            jb = jb + 1
c           
            ia = 0
            do syma = 0,7
              symi = ieor(syma,symia)
              do a = v_sym(1,syma,spina),
     $             v_sym(2,syma,spina)
                if (ia.ne.ov_off(a,symia,spini,
     $               spina))
     $               call errquit('bad ia',0, INPUT_ERR)
                do i = o_sym(1,symi,spini),
     $               o_sym(2,symi,spini)
                  ia = ia + 1
c                  write(6,*)' ia, jb, a, b, i, j, t2(a,b,i,j) ',
c     &                        ia, jb, a, b, i, j, t2(a,b,i,j)
                  if(oflip)then
                    t2buf(ia,jb) = t2(b,a,j,i)
                  else
                    t2buf(ia,jb) = t2(a,b,i,j)
                  endif
                end do
              end do
            end do
          end do
        end do
      end do
      return
      end
      subroutine pack_t1buf(t1, t1buf,
     &                      lent1buf, symjb, jb, 
     &                      olo1, ohi1, vlo1, vhi1, 
     &                      spinb, spinj)
      implicit none
c
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "amplitudes.fh"
c
      integer lent1buf
      integer symjb
      integer jb, j, b
      integer symj, symb, spinj, spinb
      integer olo1, ohi1, vlo1, vhi1
      double precision t1(vlo1:vhi1,olo1:ohi1)
      double precision t1buf(lent1buf)
      call dfill (lent1buf, 0.0d0, t1buf,1)
      jb = 0
      do symb = 0,7
        symj = ieor(symb,symjb)
        do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
          do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
            jb = jb + 1
            t1buf(jb) = t1(b,j)
          end do
        end do
      enddo
      return
      end
c $Id$
