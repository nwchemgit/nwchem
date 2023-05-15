      subroutine int_2e_sifs_b( basis, tol2e, q4, iat, jat, kat, lat, 
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi, aoint2,  ninfo, 
     $     info, ifmt, ibvtyp, values,  ibitv, buffer, clabs,
     &     ilabs, jlabs, klabs, llabs, eris )
c     
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "cfock.fh"
#include "util.fh"
#include "cscfps.fh"
#include "stdio.fh"
 
      integer ilabs(*), jlabs(*), klabs(*), llabs(*)
      double precision eris(*)
      integer basis
      integer nfock
      double precision tol2e, dentol
      double precision sijkl
      double precision q4
      integer iat, jat, kat, lat
      integer ilo, jlo, klo, llo 
      integer ihi, jhi, khi, lhi
 
      integer ncomputed, nwritten, num
      common/sif2e/ncomputed, nwritten, num
 
      logical oij, okl, oikjl, status
      integer ish, jsh, ksh, lsh, mem2, max2e
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer ibflo, jbflo, kbflo, lbflo, ibfhi, jbfhi, kbfhi, lbfhi
      integer ishbflo, ishbfhi, jshbflo, jshbfhi
      integer kshbflo, kshbfhi, lshbflo, lshbfhi
      integer jshtop, kshtop, lshtop
      integer l_scr, k_scr
 
      double precision sij, smax, denmax
      integer neri
      integer i
      integer*4 aoint2, ninfo, ifmt, ibvtyp, last, itypea, itypeb,
     &          iwait, nrec, reqnum, ierr, ibitv, nipv
      parameter(nipv=4)
      integer*4 clabs(nipv,*)
      integer*4 info(ninfo)
      double precision values(*), buffer(*)
      integer*4 msame,nmsame,nomore
      parameter(msame=0, nmsame=1, nomore= 2)
      integer n2max

c     iat, jat, kat, lat externally loop in canonical order over
c     atoms.  
c     
c     scratch space needed for 2-electron routines
c     
c     l_scr ... integral scratch space
c     l_g   ... buffer for eris
c     l_i, j, k, l ... integral labels
c     
c     Shell and function ranges for each atom
c     
cgk debug
*     write(*,*)'gk: entered int_2e_sifs_b'
cgk end
*     oerilabel=.true.
      status = .true.
      status = status .and. bas_ce2cnr(basis, iat, ishlo, ishhi)
      status = status .and. bas_ce2cnr(basis, jat, jshlo, jshhi)
      status = status .and. bas_ce2cnr(basis, kat, kshlo, kshhi)
      status = status .and. bas_ce2cnr(basis, lat, lshlo, lshhi)
      if (.not. owritefile) then
         status = status .and. bas_ce2bfr(basis, iat, ibflo, ibfhi)
         status = status .and. bas_ce2bfr(basis, jat, jbflo, jbfhi)
         status = status .and. bas_ce2bfr(basis, kat, kbflo, kbfhi)
         status = status .and. bas_ce2bfr(basis, lat, lbflo, lbfhi)
      endif
      if (.not. status) call errquit
     $     ('int_2e_sifs_b: bas_ce to cnr/bfr failed', 0, BASIS_ERR)

      call int_mem_2e4c(max2e, mem2)
      mem2 = max(mem2,max2e)

      status = .true.
      status = status .and. 
     &   ma_push_get(MT_DBL, mem2, 'int_2e_sifs_b: scr', l_scr, k_scr)
      if (.not. status) call errquit('int_2e_sifs_b: not enuf memory',0,
     &       MA_ERR)
c     
c     Permutation symmetry boolean
c     
      oij = iat .eq. jat
      okl = kat .eq. lat
      oikjl = iat.eq.kat .and. jat.eq.lat
c
c     loop over shells on each atom
c     
      smax = schwarz_max()
      n2max=info(5)
      iwait=1
      itypea=3
      itypeb=0
cgk debug
*     write(*,*)'gk: at top of loop over shell quartets'
cgk end
      do ish = ishlo, ishhi
         if (.not. bas_cn2bfr(basis, ish, ishbflo, ishbfhi))
     $        call errquit('int_2e_sifs_b: bas_cn2bfr', ish, BASIS_ERR)
         jshtop = jshhi
         if (oij) jshtop = ish
         do jsh = jshlo, jshtop
            sij = schwarz_shell(ish,jsh)
*           if (sij*smax.ge.tol2e) then
               if (.not. bas_cn2bfr(basis, jsh, jshbflo, jshbfhi))
     $              call errquit('int_2e_sifs_b: bas_cn2bfr', jsh,
     &       BASIS_ERR)
               kshtop = kshhi
               if (oikjl) kshtop = ish
               do ksh = kshlo, kshtop
                  if (.not. bas_cn2bfr(basis, ksh, kshbflo, kshbfhi))
     $                 call errquit('int_2e_sifs_b: bas_cn2bfr', ksh,
     &       BASIS_ERR)
                  lshtop = lshhi
                  if (okl) lshtop = ksh
                  if (oikjl .and. ksh.eq.ish) lshtop = jsh
                  do lsh = lshlo, lshtop
*                    if (sij*schwarz_shell(ksh,lsh)
*    $                    .ge. tol2e) then
                        if (.not. bas_cn2bfr(basis, lsh,
     $                       lshbflo, lshbfhi)) call errquit
     $                       ('int_2e_sifs_b: bas_cn2bfr', lsh,
     &       BASIS_ERR)
c     
c     Compute the shell block of integrals with labels
c     
cgk debug
*     write(*,*)'gk: calling int_2e4c'
cgk end
                        call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                       mem2, dbl_mb(k_scr), max2e, eris)
cgk debug
*     write(*,*)'gk: oerilable = ', oerilabel 
cgk end
                        call fock_2e_shell_to_labels( tol2e,
     $                       ishbflo, ishbfhi, jshbflo, jshbfhi,
     $                       kshbflo, kshbfhi, lshbflo, lshbfhi,
     &                       ilabs, jlabs, klabs, llabs, eris,
     $                       dbl_mb(k_scr),  neri)
                        ncomputed = ncomputed + neri
cgk debug
*     write(*,*)'neri = ', neri
                       do i=1,neri
cgk debug
*       write(*,2223) ilabs(i), jlabs(i), klabs(i), llabs(i), eris(i)
 2223 format(4i4,f20.12)
cgk end

                         if(num.eq.n2max) then
                           nwritten=nwritten+num
                           call sifew2( aoint2, info, nipv, num,  
     &                                  msame, itypea, itypeb, ifmt,
     &                                  ibvtyp, values, clabs, ibitv,
     &                                  buffer, iwait, nrec, reqnum,
     &                                  ierr )
                           nwritten=nwritten-num
                         endif
                         num=num+1
                         clabs(1,num)=ilabs(i)
                         clabs(2,num)=jlabs(i)
                         clabs(3,num)=klabs(i)
                         clabs(4,num)=llabs(i)
                         values(num)=eris(i)
                       enddo
*                    end if
                  end do
               end do
*           end if
         end do
      end do
c     
      status = ma_pop_stack(l_scr)
      if (.not. status) call errquit
     $     ('int_2e_sifs_b: pop stack failed?', 0, MA_ERR)
cgk debug
*     write(*,*)'gk: nwritten = ', nwritten
*     write(*,*)'gk: num = ', num
*     write(*,*)'gk: leaving int_2e_sifs_b'
cgk end
c     
      end
