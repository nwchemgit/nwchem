        subroutine mxtrn(mxin,mxout,idim,jdim)
        implicit none
        double precision mxin(idim,jdim),mxout(jdim,idim)
        integer idim,jdim
        integer i,j
c
        do i=1,idim
         do j=1,jdim
          mxout(j,i)=mxin(i,j)
         enddo
        enddo
c
        return
        end


    
      subroutine prtfmat(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,nvir)
      implicit none
      double precision f1n(nvir,nvir),f1t(nvir,nvir),f2n(nvir,nvir),
     &                 f2t(nvir,nvir),f3n(nvir,nvir),f3t(nvir,nvir),
     &                 f4n(nvir,nvir),f4t(nvir,nvir)
      integer nvir
c
      integer c,b
c
      print *,' printing out f matrices'
      do c=1,nvir
      do b=1,nvir
      write(6,100)b,c,f1n(b,c),f1t(b,c),f2n(b,c),f2t(b,c),
     &                f3n(b,c),f3t(b,c),f4n(b,c),f4t(b,c)
  100 format(2i5,4e14.5,/,10x,4e14.5)
      enddo
      enddo
c
      return
      end      
