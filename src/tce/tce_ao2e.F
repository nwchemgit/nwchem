      subroutine tce_ao2e(rtdb,d_ao2e)
c
c $Id$
c
c     Returns a SF handle for the AO 2-e integral matrix.
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "inp.fh"
      integer rtdb                 ! Run-time database
      logical nodezero             ! True if node zero
      integer d_ao2e               ! AO integrals
      integer joff(max_size)       ! j offset
      integer loff(max_size)       ! l offset
      integer work1,work2          ! Work array sizes
      integer l_work1,k_work1      ! Work array 1
      integer l_work2,k_work2      ! Work array 2
      integer l_ijkl,k_ijkl        ! AO integrals
      integer l_ij,k_ij            ! AO integrals
      integer l_pj,k_pj            ! Quarter transformed
      integer l_pq,k_pq            ! Half transformed
      integer l_pq_packed,k_pq_packed! Packed
      integer g1b,g2b
      integer sf_size_pq
      integer sf_size_kl
      integer ish,jsh,ksh,lsh
      integer i,j,k,l,kl
      integer ikjl,ijkl
      integer ilo,ihi,jlo,jhi
      integer klo,khi,llo,lhi
      integer irange,jrange
      integer krange,lrange
      integer kloff
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      integer next,count,nprocs
      integer N
      character*255 filename
      double precision dbl_size
c
c     ===========
c     Set offsets
c     ===========
c
      if (ipol.eq.1) then
        N = (nmo(1) - nfc(1) - nfv(1)) * 2
      else
        N = nmo(1) + nmo(2) - nfc(1) - nfc(2) - nfv(1) - nfv(2)
      endif
      nodezero = (ga_nodeid().eq.0)
      do j = 1, nbf
        joff(j) = (j-1)*nbf
      enddo
      do l = 1, nbf
        loff(l) = (l-1)*nbf
      enddo
c
c     ===========
c     Allocate SF
c     ===========
c
      sf_size_pq = 0
      dbl_size = 0.0d0
      do g1b = 1,noa+nob+nva+nvb
        do g2b = g1b,noa+nob+nva+nvb
          sf_size_pq = sf_size_pq + int_mb(k_range+g1b-1)
     1                            * int_mb(k_range+g2b-1)
          dbl_size   = dbl_size   + dble(int_mb(k_range+g1b-1))
     1                            * dble(int_mb(k_range+g2b-1))
        enddo
      enddo
#if !defined(EXT_INT)
      if ((dbl_size .gt. 16.0d0*1024.0d0**3).and.(ioalg.eq.0)) then
        write(LuOut,*) '--------------------------------------------'
        write(LuOut,*) 'Abort: 2 GB file size limit will be exceeded'
        write(LuOut,*) 'Use GA/SF scheme or OS large file support'
        write(LuOut,*) '--------------------------------------------'
        call errquit('tce_ao2e: abort',0,UNKNOWN_ERR)
      endif
      if (dbl_size .gt. 2.0d0*1024.0d0**3) then
        write(LuOut,*) '---------------------------------------------'
        write(LuOut,*) 'Abort: 32 bit integer overflow is anticipated'
        write(LuOut,*) 'Use a machine with a 64 bit operating system'
        write(LuOut,*) '---------------------------------------------'
        call errquit('tce_ao2e: abort',0,UNKNOWN_ERR)
      endif
#endif
      sf_size_kl = nbf**2
      if (nodezero.and.util_print('ao2e',print_default)) then
        write(LuOut,*)
        write(LuOut,9000) sf_size_pq*sf_size_kl
      endif
 9000 format(1x,'2-e (intermediate) file size = ',i15)
      call tce_filename('v2i',filename)
      if (nodezero.and.util_print('mo1e',print_default))
     1  write(LuOut,9010) filename(1:inp_strlen(filename))
 9010 format(1x,'2-e (intermediate) file name = ',A)
      call createfile(filename,d_ao2e,sf_size_pq*sf_size_kl)
c
c     =======================
c     Initialize integral API
c     =======================
c
      call int_mem_2e4c(work1,work2)
      if (.not.ma_push_get(mt_dbl,work1,'work1',l_work1,k_work1))
     1  call errquit('tce_ao2e: MA problem work1',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,work2,'work2',l_work2,k_work2))
     1  call errquit('tce_ao2e: MA problem work2',1,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf*nbf,'ij',l_ij,k_ij))
     1  call errquit('tce_ao2e: MA problem k_ij',4,MA_ERR)
      if (.not.ma_push_get(mt_dbl,N*nbf,'pj',l_pj,k_pj))
     1  call errquit('tce_ao2e: MA problem p_ij',5,MA_ERR)
      if (.not.ma_push_get(mt_dbl,N*N,'pq',l_pq,k_pq))
     1  call errquit('tce_ao2e: MA problem k_pq',6,MA_ERR)
      if (.not.ma_push_get(mt_dbl,sf_size_pq,'pq_packed',
     1  l_pq_packed,k_pq_packed)) 
     2  call errquit('tce_ao2e: MA problem pq_packed',7,
     3  MA_ERR)
c
c     ================
c     Loop over shells
c     ================
c
c      print*,'======================================'
c      print*,'nshells = ',nshells
c      print*,'======================================'
c
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do lsh=1,nshells
        do ksh=1,nshells
          if (next.eq.count) then
            if (util_print('ao2e',print_debug))
     1        write(LuOut,9190) lsh,ksh,ga_nodeid()
            if (.not.bas_cn2bfr(ao_bas_han,ksh,klo,khi))
     1        call errquit('tce_ao2e: basis fn range problem 1',0,
     2        BASIS_ERR)
            krange = khi - klo + 1
            if (.not.bas_cn2bfr(ao_bas_han,lsh,llo,lhi))
     1        call errquit('tce_ao2e: basis fn range problem 2',1,
     2        BASIS_ERR)
            lrange = lhi - llo + 1
            if (.not.ma_push_get(mt_dbl,sf_size_kl*krange*lrange,
     1        'ao2e',l_ijkl,k_ijkl))
     2        call errquit('tce_ao2e: MA problem k_ijkl',18,MA_ERR)
            do ijkl = 1,sf_size_kl*krange*lrange
              dbl_mb(k_ijkl+ijkl-1) = 0.0d0
            enddo
            do jsh=1,nshells
              do ish=1,nshells
                if (schwarz_shell(ish,ksh)*schwarz_shell(jsh,lsh) 
     1            .ge. tol2e) then
                  if (.not.bas_cn2bfr(ao_bas_han,ish,ilo,ihi))
     1              call errquit('tce_ao2e: basis fn range problem',2,
     2              BASIS_ERR)
                  irange = ihi - ilo + 1
                  if (.not.bas_cn2bfr(ao_bas_han,jsh,jlo,jhi))
     1              call errquit('tce_ao2e: basis fn range problem',3,
     2              BASIS_ERR)
                  jrange = jhi - jlo + 1
                  call int_2e4c(ao_bas_han,ish,ksh,ao_bas_han,jsh,lsh,
     1              work2,dbl_mb(k_work2),work1,dbl_mb(k_work1))
c--------------------------------------------------------------------
                  kl = 0
                  do l = llo,lhi
                    do k = klo,khi
                      kl = kl + 1
                      kloff = (kl-1)*sf_size_kl
                      do j = jlo,jhi
                        do i = ilo,ihi
                          ikjl = (((i-ilo)*krange+(k-klo))*jrange
     1                           +(j-jlo))*lrange+(l-llo)
                          dbl_mb(k_ijkl+kloff+joff(j)+i-1) = 
     1                      dbl_mb(k_work1+ikjl)
                        enddo
                      enddo
                    enddo
                  enddo
c--------------------------------------------------------------------
                endif
              enddo
            enddo
c
c     Write to SF
c
            kl = 0
            do l = llo,lhi
              do k = klo,khi
                kl = kl + 1
                kloff = (kl-1)*sf_size_kl
c
c     Half transformation
c
                call dcopy(sf_size_kl,dbl_mb(k_ijkl+kloff),1,
     1            dbl_mb(k_ij),1)
                call dgemm('T','N',N,nbf,nbf,
     1            1.0d0,dbl_mb(k_movecs_sorted),nbf,
     2            dbl_mb(k_ij),nbf,
     3            0.0d0,dbl_mb(k_pj),N)
                call dgemm('N','N',N,N,nbf,
     1            1.0d0,dbl_mb(k_pj),N,
     2            dbl_mb(k_movecs_sorted),nbf,
     3            0.0d0,dbl_mb(k_pq),N)
                call ma_pack(dbl_mb(k_pq),dbl_mb(k_pq_packed),N)
                if (util_print('ao2e',print_debug)) then
                  write(LuOut,9200) l,k
                  call ma_print_compact(dbl_mb(k_pq_packed),
     1              sf_size_pq,1,'Half transformed MO 2e')
                endif
                call put_block(d_ao2e,dbl_mb(k_pq_packed),
     1            sf_size_pq,(loff(l)+k-1)*sf_size_pq)
              enddo
            enddo
            if (.not.ma_pop_stack(l_ijkl)) 
     1        call errquit('tce_ao2e: MA problem x1',16,MA_ERR)
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call reconcilefile(d_ao2e,sf_size_pq*sf_size_kl)
c
c     ==========
c     Deallocate
c     ==========
c
      if (.not.ma_pop_stack(l_pq_packed)) 
     1  call errquit('tce_ao2e: MA problem',8,MA_ERR)
      if (.not.ma_pop_stack(l_pq)) 
     1  call errquit('tce_ao2e: MA problem',9,MA_ERR)
      if (.not.ma_pop_stack(l_pj)) 
     1  call errquit('tce_ao2e: MA problem',10,MA_ERR)
      if (.not.ma_pop_stack(l_ij)) 
     1  call errquit('tce_ao2e: MA problem',11,MA_ERR)
      if (.not.ma_pop_stack(l_work2)) 
     1  call errquit('tce_ao2e: MA problem',14,MA_ERR)
      if (.not.ma_pop_stack(l_work1)) 
     1  call errquit('tce_ao2e: MA problem',15,MA_ERR)
      if (nodezero.and.util_print('time',print_debug)) then
        write(LuOut,8000) 'createfile   ', cpusecs(1), cpusecs(51)
        write(LuOut,8000) 'deletefile   ', cpusecs(2), cpusecs(52)
        write(LuOut,8000) 'put_block    ', cpusecs(3), cpusecs(53)
        write(LuOut,8000) 'add_block    ', cpusecs(4), cpusecs(54)
        write(LuOut,8000) 'get_block    ', cpusecs(5), cpusecs(55)
        write(LuOut,8000) 'reconcilefile', cpusecs(6), cpusecs(56)
      endif
c
c     ======
c     Return
c     ======
c
 8000 format(1x,'cpu & wall time in seconds for ',A,2f20.1)
 9190 format(/,'k- and l-shells',2i3,' has been assigned to node ',i3)
 9200 format(/,'k and l',2i3)
      return
      end
