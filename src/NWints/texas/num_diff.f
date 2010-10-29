* $Id$
      program numder
      implicit real*8 (a-h,o-z)
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
c
      parameter (n=10000) 
cccc  character*16 name_p(n),name_m(n)
      character*16 name_p   ,name_m  , name_2
      character*2 direction
cccc  dimension del_plus(n,3), del_m(n,3)
c
      open(unit=7,file='delta_p',status='old',form='formatted')
      open(unit=8,file='delta_m',status='old',form='formatted')
cccc  open(unit=33,file='anali_2',status='old',form='unformatted')
c
      open(unit=34,file='compare',status='new',form='formatted')
c
      write(*,*)' enter number of points:'
      read(*,*) npoints
      write(*,*)' second differentiation over :'
      read(*,*) direction
c
      do ipoint=1,npoints
         read(7,69) name_p, icf_p,jcf_p,kcf_p,lcf_p,
     *                        der1_p_x,der1_p_y,der1_p_z
         read(8,69) name_m, icf_m,jcf_m,kcf_m,lcf_m,
     *                        der1_m_x,der1_m_y,der1_m_z
   69    format(a16,4i2,1x,3(f12.7,2x))
c
         if(name_p.eq.name_m) then
            der2_x=100.d0*(der1_p_x-der1_m_x)
            der2_y=100.d0*(der1_p_y-der1_m_y)
            der2_z=100.d0*(der1_p_z-der1_m_z)
c
            icf=icf_p
            jcf=jcf_p
            kcf=kcf_p
            lcf=lcf_p
c
            if(name_p.eq.'d /dAi   : ijkl=') then
               if(direction.eq.'Ax') name_2='d2/dAidAj: ijkl='
               if(direction.eq.'Ay') name_2='d2/dAidAj: ijkl='
               if(direction.eq.'Az') name_2='d2/dAidAj: ijkl='
c
               if(direction.eq.'Bx') name_2='d2/dAidBj: ijkl='
               if(direction.eq.'By') name_2='d2/dAidBj: ijkl='
               if(direction.eq.'Bz') name_2='d2/dAidBj: ijkl='
c
               if(direction.eq.'Cx') name_2='d2/dAidCj: ijkl='
               if(direction.eq.'Cy') name_2='d2/dAidCj: ijkl='
               if(direction.eq.'Cz') name_2='d2/dAidCj: ijkl='
c
               if(direction.eq.'Dx') name_2='d2/dAidDj: ijkl='
               if(direction.eq.'Dy') name_2='d2/dAidDj: ijkl='
               if(direction.eq.'Dz') name_2='d2/dAidDj: ijkl='
            endif
            if(name_p.eq.'d /dBi   : ijkl=') then
               if(direction.eq.'Ax') name_2='d2/dAidBj: ijkl='
               if(direction.eq.'Ay') name_2='d2/dAidBj: ijkl='
               if(direction.eq.'Az') name_2='d2/dAidBj: ijkl='
c
               if(direction.eq.'Bx') name_2='d2/dBidBj: ijkl='
               if(direction.eq.'By') name_2='d2/dBidBj: ijkl='
               if(direction.eq.'Bz') name_2='d2/dBidBj: ijkl='
c
               if(direction.eq.'Cx') name_2='d2/dBidCj: ijkl='
               if(direction.eq.'Cy') name_2='d2/dBidCj: ijkl='
               if(direction.eq.'Cz') name_2='d2/dBidCj: ijkl='
c
               if(direction.eq.'Dx') name_2='d2/dBidDj: ijkl='
               if(direction.eq.'Dy') name_2='d2/dBidDj: ijkl='
               if(direction.eq.'Dz') name_2='d2/dBidDj: ijkl='
            endif
            if(name_p.eq.'d /dCi   : ijkl=') then
               if(direction.eq.'Ax') name_2='d2/dAidCj: ijkl='
               if(direction.eq.'Ay') name_2='d2/dAidCj: ijkl='
               if(direction.eq.'Az') name_2='d2/dAidCj: ijkl='
c
               if(direction.eq.'Bx') name_2='d2/dBidCj: ijkl='
               if(direction.eq.'By') name_2='d2/dBidCj: ijkl='
               if(direction.eq.'Bz') name_2='d2/dBidCj: ijkl='
c
               if(direction.eq.'Cx') name_2='d2/dCidCj: ijkl='
               if(direction.eq.'Cy') name_2='d2/dCidCj: ijkl='
               if(direction.eq.'Cz') name_2='d2/dCidCj: ijkl='
c
               if(direction.eq.'Dx') name_2='d2/dCidDj: ijkl='
               if(direction.eq.'Dy') name_2='d2/dCidDj: ijkl='
               if(direction.eq.'Dz') name_2='d2/dCidDj: ijkl='
            endif
            if(name_p.eq.'d /dDi   : ijkl=') then
               if(direction.eq.'Ax') name_2='d2/dAidDj: ijkl='
               if(direction.eq.'Ay') name_2='d2/dAidDj: ijkl='
               if(direction.eq.'Az') name_2='d2/dAidDj: ijkl='
c
               if(direction.eq.'Bx') name_2='d2/dBidDj: ijkl='
               if(direction.eq.'By') name_2='d2/dBidDj: ijkl='
               if(direction.eq.'Bz') name_2='d2/dBidDj: ijkl='
c
               if(direction.eq.'Cx') name_2='d2/dCidDj: ijkl='
               if(direction.eq.'Cy') name_2='d2/dCidDj: ijkl='
               if(direction.eq.'Cz') name_2='d2/dCidDj: ijkl='
c
               if(direction.eq.'Dx') name_2='d2/dDidDj: ijkl='
               if(direction.eq.'Dy') name_2='d2/dDidDj: ijkl='
               if(direction.eq.'Dz') name_2='d2/dDidDj: ijkl='
            endif
c
cc             write(34,69) name_2,icf,jcf,kcf,lcf,
cc   *                      der2_x,der2_y,der2_z
               write(34,70) name_2,icf,jcf,kcf,lcf,
     *                      der2_x,der2_y,der2_z
   70    format(a16,4i2,1x,3(f12.4,2x))
c
         endif
      enddo
c
      stop
      end
