       subroutine cosmo_lebedev(iang,n,w,xyz)
       implicit none
#include "errquit.fh"
       integer iac_max, mx_l_ang
       parameter (iac_max=29, mx_l_ang=5810)

       integer iang,jang,i,n
       integer idum(mx_l_ang)
       double precision xyz(3,*),w(*)
       double precision wleb(mx_l_ang,iac_max)
       double precision zleb(3,mx_l_ang,iac_max)

       integer, external :: xc_leb_acc2nrd

       if (iang.lt.1)
     $   call errquit('cosmo_lebedev: iangleb is <1',iang,INPUT_ERR)
       if (iang.gt.iac_max)
     $   call errquit('cosmo_lebedev: iangleb is >29',iang,INPUT_ERR)

       select case (iang)
       case (3,4)
         jang = 4
       case (9,10,11)
         jang = 11
       case default
         jang = iang
       end select

       call xc_lebedev(n, .false., wleb, zleb, idum, idum, jang)
       n = xc_leb_acc2nrd(jang)

       call dcopy(3*n,zleb(1,1,jang),1,xyz,1)
       call dcopy(n,wleb(1,jang),1,w,1)

       end
