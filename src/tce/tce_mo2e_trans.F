      subroutine tce_mo2e_trans(rtdb,d_v2,kax_v2_alpha_offset,size_2e)
c
c
c
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb                 ! Run-time database
      logical nodezero             ! True if node zero
      integer d_v2                 ! MO integrals
      integer kax_v2_alpha_offset  ! MO integrals offset
      integer size_2e              ! 2e file size
c keep adding here
      integer work1,work2          ! Work array sizes
      integer l_work1,k_work1      ! Work array 1
      integer l_work2,k_work2      ! Work array 2
      integer mux,nux 
      integer mu,nu,rho,sigma
      integer mu_lo,mu_hi
      integer nu_lo,nu_hi
      integer rho_lo,rho_hi
      integer sigma_lo,sigma_hi
      integer mu_range
      integer nu_range
      integer rho_range
      integer sigma_range
      integer mu1,nu1,rho1,sigma1
      integer imu,inu,imu1,inu1
c
      integer s_aaoo,l_aaoo,k_aaoo
      integer s_aaao,l_aaao,k_aaao
c
      integer s_msns,l_msns,k_msns
      integer irho,isigma
c
      integer i,j,k,l,N
      integer k_mnrs,k_mnrj,k_mnij
      integer pair_ij,pair_kl
      integer size_k,size_l,size_kl
      integer l_klij,k_klij
      integer ipos1,ipos2
      integer iorbk,iorbl,iorbi,iorbj
      integer irow,icol,ires,addr
      integer offset_ijkl,size_ijkl,length
      integer size_i,size_j
      integer l_g1234,k_g1234
      integer l_ggmj,k_ggmj,size_ggmj
      integer size_ggij,l_ggij,k_ggij
c
      integer iha,ihb !number of corr. alpha, beta holes
      integer ipa,ipb !number of corr. alpha, beta particles
      integer g1b,g2b,g3b,g4b
c
      integer k_ind1,k_ind2
      integer INDEX_PAIR
      integer k_gpair,l_gpair,len_pair
      integer point_pair,size_g3g4
      integer l_p34,k_p34
      integer xoffset_p12,xoffset_p34
c
      integer d_v2int
      integer size_2eint
      integer len_apair
      integer l_ap_offset,k_ap_offset
      integer offset_aux,size_aux,i_aux,j_aux
      integer ixi,jxi,irshift
      integer size_mug3ij,l_mug3ij,k_mug3ij
      character*255 filename
c new offseting for 13 and 14
      integer k_aux,l_aux
c
      integer nxtask
      external nxtask
      integer next,count,nprocs
      integer Nb
      logical transpose
c --- debug ----
c      double precision xxx
c -------------- 
c 
      nodezero=(ga_nodeid().eq.0)
c 
        irshift = 0
c
c this module is called only if intorb = .true.
c N is the number of correlated orbitals
        N = nmo(1) - nfc(1) - nfv(1)
        iha = nocc(1)-nfc(1)
        ihb = nocc(ipol)-nfc(ipol)
        ipa = nmo(1)-nocc(1)-nfv(1)
        ipb = nmo(ipol)-nocc(ipol)-nfv(ipol)
c
      call tce_mo2e_offset_intorb_old(l_aux,
     &         k_aux,size_aux)
c
c     =======================
c     Initialize integral API
c     =======================
c
c ------------- HALF TRANSFORMED INTEGRALS FIRST
c ------------- ( mu >= nu | i >= j )
c
c open ga for halp transformed integrals d_v2int
c
c
c
      size_2eint=0
      do nu = 1,nshells
      do mu = nu,nshells
c remove it        i=i+1
          if (.not.bas_cn2bfr(ao_bas_han,nu,nu_lo,nu_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
          nu_range = nu_hi - nu_lo + 1
          if (.not.bas_cn2bfr(ao_bas_han,mu,mu_lo,mu_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
          mu_range = mu_hi - mu_lo + 1
       size_g3g4=mu_range*nu_range*N*(N+1)/2
       size_2eint=size_2eint+size_g3g4
      enddo
      enddo
c     size_2eint already defined
      call tce_filename('v2int',filename)
      call createfile(filename,d_v2int,size_2eint)
      call reconcilefile(d_v2int,size_2eint)
      call ga_zero(d_v2int)
c
 1122 format(1x,A,' file size   = ',i16)
      if(nodezero) 
     &write(6,1122)'v2int',size_2eint 
c
c
      len_apair=(nshells*(nshells+1))/2
      IF (.not.MA_PUSH_GET(mt_int,2*len_apair+1,'noname',
     &l_ap_offset,k_ap_offset))
     &CALL ERRQUIT('tce_mo2e_offset',0,MA_ERR)
      do i = 1,2*len_apair+1
       int_mb(k_ap_offset+i-1) = 0
      enddo
      i=0
      xoffset_p12=0
      int_mb(k_ap_offset)=len_apair
      do nu = 1,nshells
      do mu = nu,nshells
       i=i+1 
          if (.not.bas_cn2bfr(ao_bas_han,nu,nu_lo,nu_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
          nu_range = nu_hi - nu_lo + 1
          if (.not.bas_cn2bfr(ao_bas_han,mu,mu_lo,mu_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
          mu_range = mu_hi - mu_lo + 1 
       int_mb(k_ap_offset+i)=xoffset_p12
       size_g3g4=mu_range*nu_range*N*(N+1)/2
       xoffset_p12=xoffset_p12+size_g3g4
       int_mb(k_ap_offset+i+len_apair)=size_g3g4
      enddo
      enddo
c
      xoffset_p12 = 0
      size_g3g4   = 0
c
c     Pair's structure of the integral file
      call tce_mo2e_pairs(l_gpair,k_gpair,len_pair)
c
      call int_mem_2e4c(work1,work2)
      if (.not.ma_push_get(mt_dbl,work1,'work1',l_work1,k_work1))
     1  call errquit('tce_ao2e: MA problem work1',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,work2,'work2',l_work2,k_work2))
     1  call errquit('tce_ao2e: MA problem work2',1,MA_ERR)
c
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
c
      do nu = 1,nshells   !--------------
      do mu = nu,nshells   !--------------
        if (next.eq.count) then      
c test
c        write(6,900)mu,nu,ga_nodeid()
c        call util_flush(6)
c getting ranges for mu-th shell 
        if (.not.bas_cn2bfr(ao_bas_han,mu,mu_lo,mu_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
        mu_range = mu_hi - mu_lo + 1
c getting ranges for nu-th shell
        if (.not.bas_cn2bfr(ao_bas_han,nu,nu_lo,nu_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
        nu_range = nu_hi - nu_lo + 1
c opening the (mu nu | nbf j) file (AAAO) (j)
      s_aaao = mu_range*nu_range*nbf*N
      if (.not.ma_push_get(mt_dbl,s_aaao,'aaao',l_aaao,k_aaao))
     1  call errquit('tce_ao2e: MA problem aaao',0,MA_ERR)
      call dfill(s_aaao, 0.0d0, dbl_mb(k_aaao), 1)
c opening the ([mu] [nu] | rho sigma) file 
      s_msns = mu_range*nu_range*nbf*nbf
      if (.not.ma_push_get(mt_dbl,s_msns,'msns',l_msns,k_msns))
     1  call errquit('tce_ao2e: MA problem msns',0,MA_ERR)
      call dfill(s_msns, 0.0d0, dbl_mb(k_msns), 1)
c
c
CBERT doing integrals mu nu proc
c
           do rho    = 1,nshells !--------------
c
c getting ranges for rho-th shell
        if (.not.bas_cn2bfr(ao_bas_han,rho,rho_lo,rho_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
        rho_range = rho_hi - rho_lo + 1
c
           do sigma  = rho,nshells !--------------
c
c Schwarz screening !!!ACHTUNG - check this part
            if (schwarz_shell(mu,nu)*schwarz_shell(rho,sigma)
     1          .ge. tol2e) then
c getting ranges for sigma-th shell
        if (.not.bas_cn2bfr(ao_bas_han,sigma,sigma_lo,sigma_hi))
     1    call errquit('tce_ao2e: basis fn range problem 1',0,
     2    BASIS_ERR)
        sigma_range = sigma_hi - sigma_lo + 1
c getting piece of atomic 2-e integrals (mu nu | lambda sigma)
c zeroing ---
        call dfill(work1, 0.0d0, dbl_mb(k_work1), 1)
        call dfill(work2, 0.0d0, dbl_mb(k_work2), 1)
c -----------
        call int_2e4c(ao_bas_han,mu,nu,ao_bas_han,rho,sigma,
     1       work2,dbl_mb(k_work2),work1,dbl_mb(k_work1))
         j=0
         do rho1    = 1,rho_range
         do sigma1  = 1,sigma_range
         do nu1     = 1,nu_range
         do mu1     = 1,mu_range
          k_mnrs=(((mu1-1)*nu_range+nu1-1)*rho_range+rho1-1)*
     1                sigma_range+sigma1
         if(rho.ne.sigma) then 
          irho   = rho_lo + rho1 - 1
          isigma = sigma_lo + sigma1 -1
       ipos1 = (((isigma-1)*nbf+irho-1)*nu_range+nu1-1)*mu_range+mu1
       ipos2 = (((irho-1)*nbf+isigma-1)*nu_range+nu1-1)*mu_range+mu1 
          dbl_mb(k_msns+ipos1-1)=dbl_mb(k_work1+k_mnrs-1)
          dbl_mb(k_msns+ipos2-1)=dbl_mb(k_work1+k_mnrs-1)
         else 
          irho   = rho_lo + rho1 - 1
          isigma = sigma_lo + sigma1 -1
       ipos1 = (((isigma-1)*nbf+irho-1)*nu_range+nu1-1)*mu_range+mu1
          dbl_mb(k_msns+ipos1-1)=dbl_mb(k_work1+k_mnrs-1)
         end if 
c
         enddo
         enddo
         enddo
         enddo
c 
            end if !Schwartz screening
c
           enddo !sigma
         enddo !rho
CBERT atomic integrals done mu nu proc
c
c tranformation from (mu nu | nbf nbf ) => (mu nu | nbf j)
c
          do isigma=1,nbf
          do irho=1,nbf 
          do i=1,iha
           j=i
          do nu1=1,nu_range
          do mu1=1,mu_range
      k_mnrj=(((j-1)*nbf+irho-1)*nu_range+nu1-1)*mu_range+mu1
      k_mnrs=(((isigma-1)*nbf+irho-1)*nu_range+nu1-1)*mu_range+mu1
      dbl_mb(k_aaao+k_mnrj-1)=dbl_mb(k_aaao+k_mnrj-1)+
     1dbl_mb(k_msns+k_mnrs-1)*dbl_mb(k_movecs_sorted+(i-1)*nbf+isigma-1)
          enddo
          enddo
          enddo
          do i=iha+ihb+1,iha+ihb+ipa
           j=i-ihb
          do nu1=1,nu_range
          do mu1=1,mu_range
      k_mnrj=(((j-1)*nbf+irho-1)*nu_range+nu1-1)*mu_range+mu1
      k_mnrs=(((isigma-1)*nbf+irho-1)*nu_range+nu1-1)*mu_range+mu1
      dbl_mb(k_aaao+k_mnrj-1)=dbl_mb(k_aaao+k_mnrj-1)+
     1dbl_mb(k_msns+k_mnrs-1)*dbl_mb(k_movecs_sorted+(i-1)*nbf+isigma-1)
          enddo 
          enddo
          enddo
          enddo
          enddo
CBERT first transformation done mu nu proc
c closing s_msns
         if (.not.ma_pop_stack(l_msns))
     1    call errquit('tce_mo2e:s_msns MA problem',16,MA_ERR)
c
c opening the (mu nu | i j) file (AAOO) (i>=j)
      s_aaoo = mu_range*nu_range*(N*(N+1)/2)
      if (.not.ma_push_get(mt_dbl,s_aaoo,'aaoo',l_aaoo,k_aaoo))
     1  call errquit('tce_ao2e: MA problem aaoo',0,MA_ERR)
      do i = 1, s_aaoo
       dbl_mb(k_aaoo+i-1) = 0.0d0
      enddo
c transformation from ([mu] [nu] | nbf j) => ([mu] [nu] | i>=j)
      do irho=1,nbf 
c i=iha j=iha
          do i      = 1,iha
          do j      = 1,i
           pair_ij = index_pair(i,j)
          do nu1=1,nu_range
          do mu1=1,mu_range
           k_mnrj  = (((j-1)*nbf+irho-1)*nu_range+nu1-1)*
     1                 mu_range+mu1
           k_mnij  = ((pair_ij-1)*nu_range+nu1-1)*mu_range+mu1
           dbl_mb(k_aaoo+k_mnij-1)=dbl_mb(k_aaoo+k_mnij-1)+
     1     dbl_mb(k_aaao+k_mnrj-1)*
     1     dbl_mb(k_movecs_sorted+(i-1)*nbf-1+irho)
          enddo
          enddo
          enddo
          enddo
c i=ipa j=iha
          do i      = iha+ihb+1,iha+ihb+ipa
          do j      = 1,iha
           pair_ij = index_pair(i-ihb,j)
          do nu1=1,nu_range
          do mu1=1,mu_range
           k_mnrj  = (((j-1)*nbf+irho-1)*nu_range+nu1-1)*
     1                 mu_range+mu1
           k_mnij  = ((pair_ij-1)*nu_range+nu1-1)*mu_range+mu1
           dbl_mb(k_aaoo+k_mnij-1)=dbl_mb(k_aaoo+k_mnij-1)+
     1     dbl_mb(k_aaao+k_mnrj-1)*
     1     dbl_mb(k_movecs_sorted+(i-1)*nbf-1+irho)
          enddo
          enddo
          enddo
          enddo
c i=ipa j=ipa
          do i      = iha+ihb+1,iha+ihb+ipa
          do j      = iha+ihb+1,i
           pair_ij = index_pair(i-ihb,j-ihb)
          do nu1=1,nu_range
          do mu1=1,mu_range
           k_mnrj  = (((j-ihb-1)*nbf+irho-1)*nu_range+nu1-1)*
     1                 mu_range+mu1
           k_mnij  = ((pair_ij-1)*nu_range+nu1-1)*mu_range+mu1
           dbl_mb(k_aaoo+k_mnij-1)=dbl_mb(k_aaoo+k_mnij-1)+
     1     dbl_mb(k_aaao+k_mnrj-1)*
     1     dbl_mb(k_movecs_sorted+(i-1)*nbf-1+irho)
          enddo
          enddo
         enddo 
         enddo 
         enddo 
c  ga_put l_aaoo
         ixi=nshells-nu+1
         jxi=nshells-mu+1
         ipos1=(nshells*(nshells+1))/2-((ixi-1)*ixi)/2-jxi+1
         xoffset_p12 = int_mb(k_ap_offset+ipos1) 
         size_g3g4   = int_mb(k_ap_offset+len_apair+ipos1)
c CHECK POINT
         if(size_g3g4.ne.mu_range*nu_range*(N*(N+1)/2)) then
           call errquit('tce_ao2e: size_g3g4 problem',16,MA_ERR)
         end if
c
         call put_block(d_v2int,dbl_mb(k_aaoo),size_g3g4,xoffset_p12)
c 
c  closing l_aaoo
           if (.not.ma_pop_stack(l_aaoo))
     1      call errquit('tce_ao2e: laaoo MA problem',16,MA_ERR)
c  closing l_aaao
         if (.not.ma_pop_stack(l_aaao))
     1    call errquit('tce_mo2e:s_msns MA problem',16,MA_ERR)
c
          next = nxtask(nprocs,1)
        endif
        count = count + 1
      enddo !mu
      enddo !nu
      next = nxtask(-nprocs,1)
c
c
      if (.not.ma_pop_stack(l_work2))
     1  call errquit('tce_ao2e: MA problem',14,MA_ERR)
      if (.not.ma_pop_stack(l_work1))
     1  call errquit('tce_ao2e: MA problem',15,MA_ERR)
c
c first half transformation ends up here
c
      call reconcilefile(d_v2int,size_2eint)
c
c
c
c
c
c
c
c
c
c
c
c
c 
c       if(nodezero) 
c     &   call ga_summarize(.true.)
c
c
c
c
c
c
c
c
c
c
c Second part of transformation starts here
c
ccx        length = int_mb(k_v2_alpha_offset)
        length = int_mb(k_aux)
        len_pair = int_mb(k_gpair) 
        addr   = 0
        point_pair = 0
c
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
c
        do g3b = 1,noa+nva   !l
        do g4b = g3b,noa+nva !k
       if (next.eq.count) then
c test
c        write(6,901)g3b,g4b,ga_nodeid()
c        call util_flush(6)
c
         size_k  = int_mb(k_range_alpha+g4b-1)
         size_l  = int_mb(k_range_alpha+g3b-1)
         size_kl = size_k*size_l*nbf*nbf
         size_ggmj=size_k*size_l*nbf*N
c ([g4b][g3b]|mu j)
         if (.not.ma_push_get(mt_dbl,size_ggmj,'ggmj',l_ggmj,k_ggmj))
     1          call errquit('tce_ao2e: MA problem ggmj',0,MA_ERR)
         call dfill(size_ggmj, 0.0d0, dbl_mb(k_ggmj), 1)
c ([g4b][g3b]|mu nu)
         if (.not.ma_push_get(mt_dbl,size_kl,'oooo',l_klij,k_klij))
     1          call errquit('tce_ao2e: MA problem klij',0,MA_ERR)
         call dfill(size_kl, 0.0d0, dbl_mb(k_klij), 1)
c
       do mux = 1,nshells ! -------
c --- mu-loop randomization -------
         if(mux+irshift.gt.nshells) then
           mu=mux+irshift-nshells
         else
           mu=mux+irshift
         end if
c ------------------------------
         if (.not.bas_cn2bfr(ao_bas_han,mu,mu_lo,mu_hi))
     1     call errquit('tce_ao2e: basis fn range problem 1',0,
     2     BASIS_ERR)
         mu_range = mu_hi - mu_lo + 1
c
        do nux = 1,nshells  ! -------
c --- nu-loop randomization -------
         if(nux+irshift.gt.nshells) then
           nu=nux+irshift-nshells
         else
           nu=nux+irshift
         end if
c ------------------------------
c ASSUMPTION MU >= NU
         if(mu.lt.nu) go to 1888
c
         if (.not.bas_cn2bfr(ao_bas_han,nu,nu_lo,nu_hi))
     1     call errquit('tce_ao2e: basis fn range problem 1',0,
     2     BASIS_ERR)
         nu_range = nu_hi - nu_lo + 1
c open l_aaoo
         s_aaoo = mu_range*nu_range*(N*(N+1)/2)
         if (.not.ma_push_get(mt_dbl,s_aaoo,'aaoo',l_aaoo,k_aaoo))
     1     call errquit('tce_ao2e_2h: MA problem aaoo',0,MA_ERR)
         call dfill(s_aaoo, 0.0d0, dbl_mb(k_aaoo), 1)
c UNDER ASSUMPTION THAT MU >=NU
          ixi=nshells-nu+1
          jxi=nshells-mu+1
          ipos1=(nshells*(nshells+1))/2-((ixi-1)*ixi)/2-jxi+1
          xoffset_p12 = int_mb(k_ap_offset+ipos1)
          size_aux   = int_mb(k_ap_offset+len_apair+ipos1)
c --------------
c CHECK POINT
         if(size_aux.ne.mu_range*nu_range*(N*(N+1)/2)) then
           call errquit('tce_ao2e: size_aux problem',16,MA_ERR)
         end if
c
         call get_block(d_v2int,dbl_mb(k_aaoo),size_aux,xoffset_p12)
c
         do k   = 1,size_k
           iorbk=int_mb(k_offset_alpha+g4b-1)+k
         do l   = 1,size_l
           iorbl=int_mb(k_offset_alpha+g3b-1)+l 
           if(iorbk.ge.iorbl) then
            pair_kl=iorbk*(iorbk-1)/2+iorbl
           else
            pair_kl=iorbl*(iorbl-1)/2+iorbk
           end if 
         do nu1 = 1,nu_range
           inu1=nu_lo+nu1-1
         do mu1 = 1,mu_range
           imu1=mu_lo+mu1-1 
           ipos2 = ((pair_kl-1)*nu_range+nu1-1)*mu_range+mu1
          if(mu.ne.nu) then
           ipos1=(((inu1-1)*nbf+imu1-1)*size_l+l-1)*
     1                size_k+k
           dbl_mb(k_klij+ipos1-1)=dbl_mb(k_klij+ipos1-1)+
     1              dbl_mb(k_aaoo+ipos2-1)

           ipos1=(((imu1-1)*nbf+inu1-1)*size_l+l-1)*
     1                size_k+k
           dbl_mb(k_klij+ipos1-1)=dbl_mb(k_klij+ipos1-1)+
     1              dbl_mb(k_aaoo+ipos2-1)
          else
           ipos1=(((inu1-1)*nbf+imu1-1)*size_l+l-1)*
     1                size_k+k
           dbl_mb(k_klij+ipos1-1)=dbl_mb(k_klij+ipos1-1)+
     1              dbl_mb(k_aaoo+ipos2-1)
          end if
c
         enddo !mu1
         enddo !nu1
         enddo !l
         enddo !k
c close l_aaoo
           if (.not.ma_pop_stack(l_aaoo))
     1      call errquit('tce_ao2e: laaoo MA problem',16,MA_ERR)
 1888   continue
c
        enddo !  nu - loop
        enddo !  mu - loop
c transfromation from ([g4b][g3b]|mu nu) to ([g4b][g3b]| mu j) here
        do imu = 1,nbf
        do inu = 1,nbf 
        do i = 1,iha
          j=i
        do l = 1,size_l
        do k = 1,size_k
          ipos1=(((j-1)*nbf+imu-1)*size_l+l-1)*size_k+k
          ipos2=(((inu-1)*nbf+imu-1)*size_l+l-1)*size_k+k 
          dbl_mb(k_ggmj+ipos1-1)=dbl_mb(k_ggmj+ipos1-1)+
     1    dbl_mb(k_klij+ipos2-1)*
     2    dbl_mb(k_movecs_sorted+(i-1)*nbf+inu-1)
         enddo
         enddo
         enddo
c
         do i=iha+ihb+1,iha+ihb+ipa
         j=i-ihb
         do l = 1,size_l
         do k = 1,size_k
          ipos1=(((j-1)*nbf+imu-1)*size_l+l-1)*size_k+k
          ipos2=(((inu-1)*nbf+imu-1)*size_l+l-1)*size_k+k
          dbl_mb(k_ggmj+ipos1-1)=dbl_mb(k_ggmj+ipos1-1)+
     1    dbl_mb(k_klij+ipos2-1)*
     2    dbl_mb(k_movecs_sorted+(i-1)*nbf+inu-1)
        enddo
        enddo
        enddo
        enddo
        enddo  
c close l_klij
           if (.not.ma_pop_stack(l_klij))
     1      call errquit('tce_ao2e: laaoo MA problem',16,MA_ERR)
c here  l_ggmj is still open
c open ([g4b][g3b]|i>=j)
         size_ggij=size_k*size_l*(N*(N+1)/2)
         if (.not.ma_push_get(mt_dbl,size_ggij,'ggij',l_ggij,k_ggij))
     1          call errquit('tce_ao2e: MA problem ggij',0,MA_ERR)
         call dfill(size_ggij, 0.0d0, dbl_mb(k_ggij), 1)
c at this point l_ggmj AND l_ggij are open
c transfromation from ([g4b][g3b]| mu j) to ([g4b][g3b]|i>=j) here
       do imu = 1,nbf
c i=iha j=iha
       do i      = 1,iha
       do j      = 1,i
           pair_ij = index_pair(i,j)
       do l = 1,size_l
       do k = 1,size_k
           ipos2  = (((j-1)*nbf+imu-1)*size_l+l-1)*
     1                 size_k+k
           ipos1  = ((pair_ij-1)*size_l+l-1)*size_k+k
           dbl_mb(k_ggij+ipos1-1)=dbl_mb(k_ggij+ipos1-1)+
     1     dbl_mb(k_ggmj+ipos2-1)*
     1     dbl_mb(k_movecs_sorted+(i-1)*nbf-1+imu)
        enddo
        enddo
        enddo
        enddo
c i=ipa j=iha
        do i      = iha+ihb+1,iha+ihb+ipa
        do j      = 1,iha
           pair_ij = index_pair(i-ihb,j)
        do l = 1,size_l
        do k = 1,size_k
           ipos2  = (((j-1)*nbf+imu-1)*size_l+l-1)*
     1                 size_k+k
           ipos1  = ((pair_ij-1)*size_l+l-1)*size_k+k
           dbl_mb(k_ggij+ipos1-1)=dbl_mb(k_ggij+ipos1-1)+
     1     dbl_mb(k_ggmj+ipos2-1)*
     1     dbl_mb(k_movecs_sorted+(i-1)*nbf-1+imu)
         enddo
         enddo
         enddo
         enddo
c i=ipa j=ipa
         do i      = iha+ihb+1,iha+ihb+ipa
         do j      = iha+ihb+1,i
           pair_ij = index_pair(i-ihb,j-ihb)
         do l = 1,size_l
         do k = 1,size_k
           ipos2  = (((j-ihb-1)*nbf+imu-1)*size_l+l-1)*
     1                 size_k+k
           ipos1  = ((pair_ij-1)*size_l+l-1)*size_k+k
           dbl_mb(k_ggij+ipos1-1)=dbl_mb(k_ggij+ipos1-1)+
     1     dbl_mb(k_ggmj+ipos2-1)*
     1     dbl_mb(k_movecs_sorted+(i-1)*nbf-1+imu)
          enddo
          enddo
          enddo
          enddo
c
         enddo !imu
c
c on l_ggij we store ([g4b][g3b]|i>=j)
c
c   calculate point_pair index here
         ixi=noa+nva-g3b+1
         jxi=noa+nva-g4b+1
         point_pair=((noa+nva)*(noa+nva+1))/2-((ixi-1)*ixi)/2-jxi+1
c
         size_g3g4=int_mb(k_gpair+point_pair)
         if (.not.ma_push_get(mt_dbl,size_g3g4,'tilep',l_p34,k_p34))
     1          call errquit('tce_ao2e: MA problem klij',0,MA_ERR)
         call dfill(size_g3g4, 0.0d0, dbl_mb(k_p34), 1)
c
c test
c        write(6,902)g3b,g4b,ga_nodeid()
c        call util_flush(6)
c
c
c DISTRIBUTION CONSISTENT WITH THE  OFFSET FOR 2E-INT IN ORBITAL FORM
c
         xoffset_p12 = 0
c
c
         addr=int_mb(k_gpair+2*len_pair+point_pair)
c
c
         do g1b = 1,noa+nva   !j
         do g2b = g1b,noa+nva !i
c
      IF (int_mb(k_spin_alpha+g3b-1)+int_mb(k_spin_alpha+g4b-1).eq.
     &int_mb(k_spin_alpha+g1b-1)+int_mb(k_spin_alpha+g2b-1)) THEN
      IF (ieor(int_mb(k_sym_alpha+g3b-1),ieor(int_mb(k_sym_alpha+g4b-1),
     &ieor(int_mb(k_sym_alpha+g1b-1),int_mb(k_sym_alpha+g2b-1)))) .eq.
     &irrep_v) THEN
      IROW=INDEX_PAIR(g4b,g3b)
      ICOL=INDEX_PAIR(g2b,g1b)
      IF(IROW.GE.ICOL) THEN
      addr = addr + 1
c
c --- check point ---
c
      IRES=INDEX_PAIR(IROW,ICOL) 
ccx      if(ires.ne.int_mb(k_v2_alpha_offset+addr)) then
      if(ires.ne.int_mb(k_aux+addr)) then
        write(6,*)'SOMETHING WENT WRONG'
        call util_flush(6)
        stop
      end if
c -------------------
c
c read offset -> open small file -> fill it  -> ga_acc
c
ccx      offset_ijkl = int_mb(k_v2_alpha_offset+length+addr)
      offset_ijkl = int_mb(k_aux+length+addr)
      size_ijkl   = int_mb(k_range_alpha+g3b-1) *
     &int_mb(k_range_alpha+g4b-1) * int_mb(k_range_alpha+g1b-1) *
     &int_mb(k_range_alpha+g2b-1)
      size_i = int_mb(k_range_alpha+g2b-1)
      size_j = int_mb(k_range_alpha+g1b-1)
c
c
c (k l | i j) (k,l,i,j - orbital indices corresponding to alpha orbitals) 
c
         do k   = 1,size_k
         do l   = 1,size_l      
         do i   = 1,size_i
         do j   = 1,size_j
c redefined definition of ipos1
          ipos1 = (((j-1)*size_i+i-1)*size_l+l-1)*size_k+k+
     1            xoffset_p12 
          iorbi = int_mb(k_offset_alpha+g2b-1)+i
          iorbj = int_mb(k_offset_alpha+g1b-1)+j          
          if(iorbi.ge.iorbj) then
           pair_ij = index_pair(iorbi,iorbj)
          else 
           pair_ij = index_pair(iorbj,iorbi)
          end if
          ipos2 = ((pair_ij-1)*size_l+l-1)*size_k+k
c
          dbl_mb(k_p34+ipos1-1)=dbl_mb(k_p34+ipos1-1)+
     1                          dbl_mb(k_ggij+ipos2-1)
         enddo
         enddo
         enddo
         enddo
c
         xoffset_p12 = xoffset_p12 + size_ijkl
c
         endif
         endif
         endif
         enddo !g2b i
         enddo !g1b j
c test
c        write(6,903)g3b,g4b,ga_nodeid()
c        call util_flush(6)
c
         xoffset_p34 = int_mb(k_gpair+len_pair+point_pair)
c
c         call ga_acc(d_v2,xoffset_p34+1,xoffset_p34+size_g3g4,1,1,
c     1               dbl_mb(k_p34),1,1.0d0) 
c
         call put_block(d_v2,dbl_mb(k_p34),size_g3g4,xoffset_p34)
c
         if (.not.ma_pop_stack(l_p34))
     1    call errquit('tce_mo2e_t: l_p34 MA problem',16,MA_ERR)
c
         if (.not.ma_pop_stack(l_ggij))
     1    call errquit('tce_mo2e_t: ggij MA problem',16,MA_ERR)
c
         if (.not.ma_pop_stack(l_ggmj))
     1    call errquit('tce_mo2e_t: ggmj MA problem',16,MA_ERR)
c
          next = nxtask(nprocs,1)
        endif
        count = count + 1
c
        enddo  !g4b k
        enddo  !g3b l
c
      next = nxtask(-nprocs,1)
c
      call reconcilefile(d_v2,size_2e)
c
c
c     ==============
c     Deallocate API
c     ==============
c
      if (.not.ma_pop_stack(l_gpair))
     1  call errquit('tce_ao2e: MA problem',15,MA_ERR)
c
      if (.not.ma_pop_stack(l_ap_offset))
     1  call errquit('tce_ao2e: MA problem',15,MA_ERR)
c
      if (.not.ma_pop_stack(l_aux))
     1  call errquit('tce_ao2e: MA problem',15,MA_ERR)
c closing d_v2int ga
      call deletefile(d_v2int) 
c
  900 format('first half mu nu CPU   ',3i7)
  901 format('second half g3b g4b CPU',3i7)
  902 format('second half after all mu-nu: g3b g4b CPU',3i7)
  903 format('second half after g1-g2:     g3b g4b CPU',3i7)
c
      return
      end
c $Id$
