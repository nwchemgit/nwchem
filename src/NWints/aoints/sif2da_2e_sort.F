      subroutine sif2da_2e_sort(aodens, l2rec, n2max, info, ibas, nbft,
     &  nshell, b2s, sminb, smaxb, frc_2el, nat, lbuf, buf, lscr,
     &  scr)

#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"

      integer*4 aodens
      integer*4 l2rec
      integer*4 n2max
      integer*4 info(*)
      integer ibas
      integer nbft
      integer nshell
      integer b2s(nbft) ! basis number to shell mapping array (output)
      integer sminb(nshell) ! min basis in shell (output)
      integer smaxb(nshell) ! max basis in shell (output)
      double precision frc_2el(3,nat)
      integer nat
      integer lbuf
      double precision buf(lbuf)
      integer lscr
      double precision scr(lscr)

c     local
      integer ish, jsh, ksh, lsh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer i, j, k, l, ij, kl
      integer ijkl, mx_d2, sh_ijkl, ijsh, klsh, sh_ijkl_len
      integer ij_max, kl_min, mxp
      integer sizeof_int, szlabs
      integer l_sq, l_d2, l_vb, l_lb, l_buf, l_dbg
      integer k_sq, k_d2, k_vb, k_lb, k_buf, k_dbg
      integer cijkl
      logical status
      integer nipv
      parameter(nipv=4)
      integer s, t, u, v, st, uv, stuv, shoff, sqindx
      integer nsp, is, js, ks, ls, a, b, c, d, ab, cd, abcd
      integer imin, jmin, kmin, lmin, imax, jmax, kmax, lmax
      integer slen, tlen, ulen, vlen, stlen
      integer i8n2max, i8l2rec
      integer avail
      integer iln, kln

cgk debug
*     write(*,*)'gk: entered sif2da_2e_sort'
cgk end


c     maximimum number of unique density matrix elements
      mx_d2 = (nbft*(nbft+1)+2)*nbft*(nbft+1)/8

c     maximum pairs
      mxp = (nbft*(nbft+1))/2
      nsp=((nshell+1)*nshell)/2
cgk debug
*     write(*,*)'gk: nsp = ', nsp
cgk end

c     generate local data to avoid overhead
      do lsh=1,nshell
       status=bas_cn2bfr(ibas,lsh,llo,lhi)
       do l=llo,lhi
        b2s(l)=lsh
       enddo
       sminb(lsh)=llo
       smaxb(lsh)=lhi  
      enddo

cgk debug
*     write(*,*)'gk: nbft = ', nbft
*     write(*,*)'gk: b2s:'
*     do l=1,nbft
*      write(*,*)l,  b2s(l)
*     enddo
*     write(*,*)'gk: sminb, smaxb:'
*     do lsh=1,nshell
*       write(*,*)lsh, sminb(lsh), smaxb(lsh)
*     enddo
cgk end

c     Get memory for shell block off set.  Diminsioned for 
c     number of shell pairs by number of shell pairs
      if (.not. ma_push_get(mt_int,nsp*nsp,'num_shell_quart',
     $     l_sq, k_sq)) call errquit
     &      ('sif2da_2e_sort;could not allocate l_sq',nsp*nsp, 
     &       MA_ERR)
cgk debug
c     Get memory for debuggin array
      if (.not. ma_push_get(mt_int,4*mx_d2,'debug_array',
     $     l_dbg, k_dbg)) call errquit
     &      ('sif2da_2e_sort;could not allocate dbg array',nsp*nsp, 
     &       MA_ERR)
cgk end

      call ifill(4*mx_d2,0,int_mb(k_dbg),1)
      call ifill(nsp*nsp, 0, int_mb(k_sq), 1)

c     go order similar to aodens and figure shell quartet offests
      ijkl = 0
      sh_ijkl=0
      sh_ijkl_len=0
      shoff=0
      do ksh=1,nshell
       klo=sminb(ksh)
       khi=smaxb(ksh)
       kln=khi-klo+1
       do lsh=1,ksh
*       klsh=nshell*(lsh-1)-((lsh-1)*lsh)/2+ksh
        klsh=((ksh-1)*ksh)/2+lsh
        llo=sminb(lsh)
        lhi=smaxb(lsh)
        do ish=ksh,nshell
         ilo=sminb(ish)
         ihi=smaxb(ish)
         iln=ihi-ilo+1
         do jsh=1,ish
*         ijsh=nshell*(jsh-1)-((jsh-1)*jsh)/2+ish
          ijsh=((ish-1)*ish)/2+jsh
          jlo=sminb(jsh)
          jhi=smaxb(jsh)
          if(ijsh.ge.klsh) then
           sh_ijkl=sh_ijkl+1
           sqindx=(klsh-1)*nsp+ijsh
           int_mb(k_sq-1+sqindx)=shoff
cgk debug
*          write(6,2323) '(',ish,',',jsh,';',ksh,',',lsh,')'
*          write(6,'(2(a,i4))')'st=',ijsh, ' uv=', klsh
*          write(6,'(2(a,i5))') 'sqindx = ', sqindx, 
*    &       ' shoff=', int_mb(k_sq-1+sqindx)
cgk end
           sh_ijkl_len=0
           do k=klo,khi
            if(lsh.eq.ksh) lhi=k
            do l=llo,lhi
             kl=((k-1)*k)/2+l
             do i=ilo,ihi
              if(jsh.eq.ish) jhi=i
              do j=jlo,jhi
               ij=((i-1)*i)/2+j
               if(ijsh.eq.klsh) then
                if(ij.ge.kl) then
                 sh_ijkl_len=sh_ijkl_len+1
                 ijkl=ijkl+1
cgk debug
                 int_mb(k_dbg-1+4*(ijkl-1)+1)=i
                 int_mb(k_dbg-1+4*(ijkl-1)+2)=j
                 int_mb(k_dbg-1+4*(ijkl-1)+3)=k
                 int_mb(k_dbg-1+4*(ijkl-1)+4)=l
cgk end
cgk debug
*                write(6,'(a,i4)') '(ijkl)=', ijkl
*                write(6,2323) '(',ish,',',jsh,';',ksh,',',lsh,')'
*                write(6,2323) '[',i,',',j,';',k,',',l,']'
cgk end
                endif
               else
                sh_ijkl_len=sh_ijkl_len+1
                ijkl=ijkl+1
cgk debug
                int_mb(k_dbg-1+4*(ijkl-1)+1)=i
                int_mb(k_dbg-1+4*(ijkl-1)+2)=j
                int_mb(k_dbg-1+4*(ijkl-1)+3)=k
                int_mb(k_dbg-1+4*(ijkl-1)+4)=l
cgk end
cgk debug
 2322 format(2(a,i2))
*               write(6,'(a,i4)') '(ijkl)=', ijkl
*               write(6,2323) '(',ish,',',jsh,';',ksh,',',lsh,')'
*               write(6,2323) '[',i,',',j,';',k,',',l,']'
 2323 format(4(a1,i2),a1)
               endif
              enddo
             enddo
            enddo
           enddo
           shoff=shoff+sh_ijkl_len
*          write(*,*)'gk: shell length=',sh_ijkl_len
          endif
         enddo
        enddo
       enddo
      enddo

cgk debug
*      write(*,*)'gk: ijkl = ', ijkl
*      write(*,*)'gk: mx_d2 = ', mx_d2

*     call MA_summarize_allocated_blocks
*     avail = MA_inquire_avail(mt_dbl)
*     write(*,*)'gk: available mem in sif2da before d2 is allocated = ',
*    &  avail
cgk end


c     memory for the *entire* 2 electron density array
      if (.not. ma_push_get(mt_dbl, mx_d2,'d2 array',l_d2,
     $     k_d2)) call errquit
     &      ('sif2da_2e_sort;could not allocate l_d2',mx_d2, MA_ERR)
c     memory for the 2 electron sifs buffer
      i8n2max=n2max
      if (.not. ma_push_get(mt_dbl,i8n2max,'vals buf array',l_vb,
     $     k_vb)) call errquit
     &      ('sif2da_2e_sort;could not allocate l_vb',i8n2max, MA_ERR)

c     get memory for sifs 2e labels
      sizeof_int=MA_sizeof(MT_INT, 1, MT_BYTE)
      if(sizeof_int.eq.4) then
        szlabs=(i8n2max+1)/2
      elseif(sizeof_int.eq.8) then
        szlabs=i8n2max
      else    
        call errquit('wrt_dft_aoints: unsupporte integer type size=',
     &    sizeof_int, MA_ERR) 
      endif   
      if (.not. MA_push_get(MT_Int,nipv*szlabs,'wrt_dft_aoints:dlab',
     &  l_lb, k_lb))
     $  call errquit('wrt_dft_aoints: ma failed for sif labs', 
     &    nipv*szlabs, MA_ERR) 
c     memory for the 2 electron sifs buffer
      i8l2rec=l2rec
      if (.not. ma_push_get(mt_dbl,i8l2rec,'vals buf array',l_buf,
     $     k_buf)) call errquit
     &      ('sif2da_2e_sort;could not allocate l_buf',l2rec, MA_ERR)

c
c     generate array of integrals from sifs file
c
cgk debug
*     write(*,*)'gk: calling sif2arr'
*     open(unit=26,form='formatted')
cgk end
      call sif2arr(aodens, dbl_mb(k_d2), mx_d2, dbl_mb(k_vb), nipv,
     & int_mb(k_lb), dbl_mb(k_buf), info, b2s, sminb, smaxb, 
     & int_mb(k_sq), nsp, nshell, int_mb(k_dbg))
cgk debug
*     write(*,*)'gk: back from sif2arr'
*     close(26)
cgk end

c     dump sif 2e buffers
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('sif2da_2e: pop failed at l_buf', 0, GA_ERR)
      if (.not. MA_pop_stack(l_lb)) call errquit
     $     ('sif2da_2e: pop failed at l_lb', 0, GA_ERR)
      if (.not. MA_pop_stack(l_vb)) call errquit
     $     ('sif2da_2e: pop failed at l_vb', 0, GA_ERR)

c     get memory for the gradient integrals
cgk debug
*     write(*,*)'gk: in sif2da, lbuf = ', lbuf
*     write(*,*)'gk: in sif2da, lscr = ', lscr
cgk end

      call d2geri_trace(mx_d2, dbl_mb(k_d2), lscr, scr, lbuf, buf,
     &  nbft, b2s, nshell, sminb, smaxb, nsp, int_mb(k_sq), nat, 
     &  frc_2el, ibas, int_mb(k_dbg))
 
cgk debug
*     write(*,*)'gk: in sif2da_2e_sort back for d2geri_trace'
cgk end


      if (.not. MA_pop_stack(l_d2)) call errquit
     $     ('sif2da_2e: pop failed at l_d2', 0, GA_ERR)
cgk debug
      if (.not. MA_pop_stack(l_dbg)) call errquit
     $     ('sif2da_2e: pop failed at l_dbg', 0, GA_ERR)
cgk end
      if (.not. MA_pop_stack(l_sq)) call errquit
     $     ('sif2da_2e: pop failed at l_sq', 0, GA_ERR)
      return 
      end

