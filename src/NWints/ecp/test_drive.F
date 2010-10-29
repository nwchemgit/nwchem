* $Id$
      implicit none
      Integer LUIn, LUOut, LUErr
      Parameter (LuIn = 5, LuOut = 6, LuErr = 0)
C###      Common /StdIO_Units/ LuIn, LuOut, LuErr
#include "ecp_consts.fh"
      integer i,j,k,l,m,n,ii,jj,lla,llb,
     &    iz_bas,ic_bas,iz_ecp,ic_ecp,jz_bas,jc_bas,n_bas,nprim,ncont,
     &    l_bas,l_ecp,i_C,n_C,idum,
     &    n_zeta_C,n_coef_C,n_C,l_ecp_max,lmax_both,
     &    n_prim_A,n_cont_A,l_A,i_c_A,nn_A,nb_A,
     &    n_prim_B,n_cont_B,l_B,i_c_B,nn_B,nb_B,
     &    n_int,n_blk,sphcart,lscr,ibug
      integer n_prim_C(0:2,-1:5,20),n_colc_C(-1:5,20),
     &    ind_z(-1:5,20),ind_c(-1:5,20),
     &    l_C(10),i_cent_C(10),n_bas_C(10),
     &    l_bas_C(20,10),nprim_bas(20,10),ncont_bas(20,10)

      integer lcsco,lscr
      parameter (lcsco = 10000,lscr=1000000)
      logical DryRun
      double precision w3,wint,wn,zz,
     &    xyz(3,10),z_bas(1000),c_bas(10000),
     &    z_ecp(1000),c_ecp(1000),
     &    scr(lscr),ecp_ints(10000)
      double precision csco(lcsco)
      character*2 type(10)
      character*1 llabel(0:2)
      data type/'  ','x ','y ','z ','x2','xy','xz','y2','yz','z2'/
      data llabel/'s','p','d'/
*
*  Read data
*
      call ecp_init_con
      w3 = one/sqrt(3.0d0)
      read (LuIn,*) ibug
      iz_bas = 0
      ic_bas = 0
      iz_ecp = 1
      ic_ecp = 1
      read (LUIn,*) n_C
      l_bas = 0
      l_ecp = 0
      do i_C = 1,n_C
        read (LUIn,*) (xyz(i,i_C),i=1,3)
        read (LUIn,*) n_bas
        n_bas_C(i_C) = n_bas
        do n = 1,n_bas
          read (LUIn,*) l,nprim,ncont
          l_bas_C(n,i_C) = l
          nprim_bas(n,i_C) = nprim
          ncont_bas(n,i_C) = ncont
          l_bas = max(l_bas,l_bas_C(n,i_C))
          do m = 1,nprim
            iz_bas = iz_bas+1
            read (LUIn,*) z_bas(iz_bas),(c_bas(ic_bas+k),k=1,ncont)
*
*         normalize basis functions
*
            zz =  z_bas(iz_bas)
            zz = zz+zz
            wn = sqrt((zz/pi)**3)
            if (l .gt. 0) wn = wn*(zz+zz)**l
            wn = sqrt(wn)
            do k = 1,ncont
              c_bas(ic_bas+k) = c_bas(ic_bas+k)*wn
            end do
            ic_bas = ic_bas+ncont
          end do
        end do
*
*     Read ECP data
*
        read (LUIn,*) l_C(i_C) 
        l_ecp = max(l_ecp,l_C(i_C))
        i_cent_C(i_C) = i_C
*
*     Local potential
*
        read (LUIn,*)  (n_prim_C(k,-1,i_C),k=0,2)
        ind_z(-1,i_C) = iz_ecp
        ind_c(-1,i_C) = ic_ecp
        n = 0
        do k = 0,2
          n = n+n_prim_C(k,-1,i_C)
          do m = 1,n_prim_C(k,-1,i_C)
            read (LUIn,*) j,z_ecp(iz_ecp),c_ecp(ic_ecp)
            if (j .ne. k) then
              write (LUOut,*) 'ECP data mismatch'
              stop
            end if
            iz_ecp = iz_ecp+1
            ic_ecp = ic_ecp+1
          end do
        end do
        n_colc_C(-1,i_C) = n
*
*     Scalar nonlocal potential
*
        do l = 0,l_C(i_C)
          read (LUIn,*)  (n_prim_C(k,l,i_C),k=0,2)
          ind_z(l,i_C) = iz_ecp
          ind_c(l,i_C) = ic_ecp
          n = 0
          do k = 0,2
            n = n+n_prim_C(k,l,i_C)
            do m = 1,n_prim_C(k,l,i_C)
              read (LUIn,*) j,z_ecp(iz_ecp),c_ecp(ic_ecp)
              if (j .ne. k) then
                write (LUOut,*) 'ECP data mismatch'
                stop
              end if
              iz_ecp = iz_ecp+1
              ic_ecp = ic_ecp+1
            end do
          end do
          n_colc_C(l,i_C) = n
        end do
*
*     Spin-orbit nonlocal potential
*
        do l = 1,l_C(i_C)
          read (LUIn,*)  (n_prim_C(k,l,i_C+n_C),k=0,2)
          ind_z(l,i_C+n_C) = iz_ecp
          ind_c(l,i_C+n_C) = ic_ecp
          n = 0
          do k = 0,2
            n = n+n_prim_C(k,l,i_C+n_C)
            do m = 1,n_prim_C(k,l,i_C+n_C)
              read (LUIn,*) j,z_ecp(iz_ecp),c_ecp(ic_ecp)
              if (j .ne. k) then
                write (LUOut,*) 'ECP data mismatch'
                stop
              end if
              iz_ecp = iz_ecp+1
              ic_ecp = ic_ecp+1
            end do
          end do
          n_colc_C(l,i_C+n_C) = n
        end do
      end do
      n_zeta_C = iz_ecp-1
      n_coef_C = ic_ecp-1
*
*  Initialize spherical transformation
*
      call defNxyz(l_bas)
      lmax_both = l_ecp + l_bas + 2 
      call ecp_init_c2s(lmax_both,csco,csco,lcsco,1,1,.true.,idum)
      write (LuOut,*) 'Cstrans needs',idum
      call ecp_init_c2s(lmax_both,csco,csco,lcsco,1,1,.false.,idum)
*
*  Call integral routine
*
      l_ecp_max = 5
      sphcart = 0
      DryRun = .false.
*
      n_blk = 1
      iz_bas = 1
      ic_bas = 1
      do i = 1, n_C
        do k = 1,n_bas_C(i)
          l_A = l_bas_C(k,i)
          lla = l_A*l_A
          n_prim_A = nprim_bas(k,i)
          n_cont_A = ncont_bas(k,i)
          nn_A = (l_A+1)*(l_A+2)/2
          nb_A = nn_A*n_cont_A
          i_c_A = i
*
          jz_bas = 1
          jc_bas = 1
          do j = 1,i
            n_bas = n_bas_C(j)
            if (j .eq. i) n_bas = k
            do m = 1,n_bas
              l_B = l_bas_C(m,j)
              llb = l_B*l_B
              n_prim_B = nprim_bas(m,j)
              n_cont_B = ncont_bas(m,j)
              nn_B = (l_B+1)*(l_B+2)/2
              nb_B = nn_B*n_cont_B
              i_c_B = j
              n_int = nb_A*nb_B
              call ecp_integral (
     &            xyz(1,i),z_bas(iz_bas),c_bas(ic_bas),
     &            n_prim_A,n_cont_A,l_A,i_c_A,
     &            xyz(1,j),z_bas(jz_bas),c_bas(jc_bas),
     &            n_prim_B,n_cont_B,l_B,i_c_B,
     &            xyz,z_ecp,c_ecp,n_prim_C,n_colc_C,
     &            ind_z,ind_c,n_zeta_C,n_coef_C,
     &            l_C,i_cent_C,n_C,l_ecp_max,
     &            sphcart,csco,lcsco,
     &            ecp_ints,n_int,n_blk,
     &            DryRun,scr,lscr,ibug)
              n = 1
              do ii = 1,nn_A
                do jj = 1,nn_B
                  wint = ecp_ints(n)
                  if (l_B .eq. 2 .and. (jj .eq. 1 .or. jj .eq. 4 .or.
     &                jj .eq. 6)) wint = wint*w3
                  if (l_A .eq. 2 .and. (ii .eq. 1 .or. ii .eq. 4 .or.
     &                ii .eq. 6)) wint = wint*w3
                  if (abs(wint) .gt. 1.0d-10)
     &                write (LuOut,600) i,llabel(l_A),k,type(ii+lla),
     &                j,llabel(l_B),m,type(jj+llb),wint
  600             format(2(i1,a1,i1,1x,a2,4x),F26.16)
                  n = n+1
                end do
              end do
*
              jz_bas = jz_bas+n_prim_B
              jc_bas = jc_bas+n_prim_B*n_cont_B
            end do
          end do
          iz_bas = iz_bas+n_prim_A
          ic_bas = ic_bas+n_prim_A*n_cont_A
        end do
      end do
*
*     Spin-orbit integrals
*
      n_blk = 3
      iz_bas = 1
      ic_bas = 1
      do i = 1, n_C
        do k = 1,n_bas_C(i)
          l_A = l_bas_C(k,i)
          lla = l_A*l_A
          n_prim_A = nprim_bas(k,i)
          n_cont_A = ncont_bas(k,i)
          nn_A = (l_A+1)*(l_A+2)/2
          nb_A = nn_A*n_cont_A
          i_c_A = i
*
          jz_bas = 1
          jc_bas = 1
          do j = 1,i
            n_bas = n_bas_C(j)
            if (j .eq. i) n_bas = k
            do m = 1,n_bas
              l_B = l_bas_C(m,j)
              llb = l_B*l_B
              n_prim_B = nprim_bas(m,j)
              n_cont_B = ncont_bas(m,j)
              nn_B = (l_B+1)*(l_B+2)/2
              nb_B = nn_B*n_cont_B
              i_c_B = j
              n_int = nb_A*nb_B
              write (6,*) '**********', i,k,j,m, '**********'
              call ecp_integral (
     &            xyz(1,i),z_bas(iz_bas),c_bas(ic_bas),
     &            n_prim_A,n_cont_A,l_A,i_c_A,
     &            xyz(1,j),z_bas(jz_bas),c_bas(jc_bas),
     &            n_prim_B,n_cont_B,l_B,i_c_B,
     &            xyz,z_ecp,c_ecp,n_prim_C,n_colc_C,
     &            ind_z,ind_c,n_zeta_C,n_coef_C,
     &            l_C,i_cent_C,n_C,l_ecp_max,
     &            sphcart,csco,lcsco,
     &            ecp_ints,n_int,n_blk,
     &            DryRun,scr,lscr,ibug)
              n = 1
              do l = 1,3
                write (LuOut,*) type(5-l),'spin-orbit'
                do ii = 1,nn_A
                  do jj = 1,nn_B
                    wint = ecp_ints(n)
                    if (l_B .eq. 2 .and. (jj .eq. 1 .or. jj .eq. 4 .or.
     &                  jj .eq. 6)) wint = wint*w3
                    if (l_A .eq. 2 .and. (ii .eq. 1 .or. ii .eq. 4 .or.
     &                  ii .eq. 6)) wint = wint*w3
                    if (abs(wint) .gt. 1.0d-10)
     &                  write (LuOut,600) i,llabel(l_A),k,type(ii+lla),
     &                  j,llabel(l_B),m,type(jj+llb),wint
                    n = n+1
                  end do
                end do
              end do
*
              jz_bas = jz_bas+n_prim_B
              jc_bas = jc_bas+n_prim_B*n_cont_B
            end do
          end do
          iz_bas = iz_bas+n_prim_A
          ic_bas = ic_bas+n_prim_A*n_cont_A
        end do
      end do
*
      stop
      end
