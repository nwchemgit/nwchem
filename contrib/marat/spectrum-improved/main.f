       program main
       implicit none
       logical ofile
       integer nmax
       parameter(nmax=20)
       double precision lmax(nmax),sigma(nmax),sum(nmax)
       double precision f(nmax),e(nmax)
       character*(180) buffer
       integer i,n
       double precision l1,l2,tmp
c
       inquire(file="spectrum.dat",exist=ofile)
       if(ofile) then
          open(10,file="spectrum.dat",
     $            form='formatted',status='old')
       else
         write(6,*) "Please provide spectrum.dat file in the format"
         return
       end if 
       n = 0
10     continue
       read(10,'(A180)',ERR=30,END=30) buffer
       n = n+1
       if(n.gt.nmax) then
         write(6,*) "increase nmax"
         return
       end if
       read(buffer,*) lmax(n),f(n),sigma(n)
       goto 10
30     close(10)
       l1=lmax(1)
       l2=lmax(1)
       do i=1,n
         write(*,*) i,lmax(i),f(i),sigma(i)
         call integrate(lmax(i),sigma(i),sum(i))
         write(*,*) "integral equals to",sum(i)
         e(i)=f(i)/(0.043*sum(i))
         write(*,*) "e=",e(i)
         tmp = lmax(i)-5*sigma(i)
         if(tmp.lt.l1) l1=tmp
         tmp = lmax(i)+5*sigma(i)
         if(tmp.gt.l2) l2=tmp
       end do

       write(*,*) "proposed boundaries",l1,l2
       call plot(n,l1,l2,lmax,sigma,e)
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

       subroutine plot(n,l1,l2,lmax,sigma,e)
       implicit none
       integer n
       double precision l1,l2
       double precision lmax(n),sigma(n),e(n)
c       
       integer i
       integer k,kmax
       double precision f,tmp
       double precision x,dx
       double precision dl

       dx = 2.0d0
       kmax = int((l2-l1)/dx)
       write(*,*) "number of points",kmax
       x = l1
       do k=1,kmax-1
         x=x+dx
         f = 0.0d0
         do i=1,n
           f =f+e(i)*exp(-((lmax(i)-x)/sigma(i))**2)
         end do
         write(33,*) x,f
       end do
       return
       
       end
       subroutine plot1(lmax,sigma,e)
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


