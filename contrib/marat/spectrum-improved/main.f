       program main
       implicit none
       double precision lmax,sigma,sum
       double precision f,e
       logical ofile
c
       inquire(file="spectrum.dat",exist=ofile)
       if(ofile) then
          open(10,file="spectrum.dat",
     $            form='formatted',status='old')
c         lattice vectors a,b,c
          read(10,*) lmax
          read(10,*) sigma
          read(10,*) f
          close(10)
       else
       write(6,*) "Please provide spectrum.dat file in the format"
       end if 
       call integrate(lmax,sigma,sum)
       write(*,*) "integral equals to",sum
       e=f/(0.043*sum)
       write(*,*) "e=",e
       call plot(lmax,sigma,e)
       end
c
       subroutine integrate(lmax,sigma,sum)
       implicit none
       double precision lmax,sigma,sum
c       
       integer i,n
       double precision f,tmp
       double precision x,dx
       double precision dl

       dx = 0.1d0
       dl=5.0d0*sigma
       n = int(2*dl/dx)
       write(*,*) "number of points",n
       x = lmax-dl
       f = exp(-((lmax-x)/sigma)**2)/x**2
       sum = 0.5d0*f
       do i=1,n-1
         x=x+dx
         f = exp(-((lmax-x)/sigma)**2)/x**2
         sum = sum + f
       end do
       x=x+dx
       f = exp(-((lmax-x)/sigma)**2)/x**2
       sum = sum +0.5d0*f
       sum = sum*2*dl/dble(n)
       return
       
       end

       subroutine plot(lmax,sigma,e)
       implicit none
       double precision lmax,sigma,e
c       
       integer i,n
       double precision f,tmp
       double precision x,dx
       double precision dl

       dx = 0.1d0
       dl=5.0d0*sigma
       n = int(2*dl/dx)
       write(*,*) "number of points",n
       x = lmax-dl
       do i=1,n-1
         x=x+dx
         f = e*exp(-((lmax-x)/sigma)**2)
         write(33,*) x,f
       end do
       return
       
       end


