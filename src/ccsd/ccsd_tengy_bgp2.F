c!
c! 1) eorb_m( 0 ) => eorb( ncor + nocc ), allos to have 1 pointer instead of 3
c! 2) #include "ccsdps.fh" is unnecessary
c! 3) constants to variables is bad 
c! 4) remove if (occsdps) for a moment
c! 
c! Original 60794725
c! Computed 42410408, intrinsics for d0-d9
        subroutine ccsd_tengy_bgp2(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     1             dintc,dintx,t1v,eorb_r,eaijk,emp4,emp5,nvir)
        implicit none
        integer ncor,nocc,nvir
        integer b,c
        double precision eaijk,emp4,emp5,denom
        double precision f1n(nvir,nvir),f1t(nvir,nvir)
        double precision f2n(nvir,nvir),f2t(nvir,nvir)
        double precision f3n(nvir,nvir),f3t(nvir,nvir)
        double precision f4n(nvir,nvir),f4t(nvir,nvir)
        double precision dintc(nvir),dintx(nvir),t1v(nvir),eorb_r(*)

        complex(8), parameter :: two = (2d0,2d0), three = (3d0,3d0)
        complex(8) af1t,af1n,af2t,af2n,af3t,af3n,af4t,af4n,
     1             adix,adiv,at1v,axmp,adic,b1,b2,b3,b4,b5,b6,b0,b9
        complex(8) c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,ca,
     1             bf1t,bf2t,bf3t,bf4t,bf1n,bf2n,bf3n
        complex(8) d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,da,
     1             cf1t,cf2t,cf3t,cf4t,cf1n,cf2n,cf3n,cf4n
        complex(8) e0,e1,e3,e5,e9,ef1t,ef2t,ef3t,ef4t,ef1n,ef2n,ef3n,
     1             bxmp,bdiv
        real(8) f0,f1

        real(8) ediv(nvir,nvir)

c! ediv is symmetric => ediv( b,c ) = ediv( c,b )
        do b = 1, nvir
           do c = b, nvir
              ediv( c,b ) = 1d0 / ( eaijk - eorb_r( b ) - eorb_r( c ) )
              ediv( b,c ) = ediv( c,b )
           enddo    
        enddo

        axmp = ( emp5,0d0 )

        do b = 1, nvir
            do c = 1, nvir, 2
                af1n = LOADFP( f1n( c,b ) )
                af1t = LOADFP( f1t( c,b ) )
                af2n = LOADFP( f2n( c,b ) )
                af2t = LOADFP( f2t( c,b ) )
                af3n = LOADFP( f3n( c,b ) )
                af3t = LOADFP( f3t( c,b ) )
                af4n = LOADFP( f4n( c,b ) )
                af4t = LOADFP( f4t( c,b ) )
                at1v = LOADFP( t1v( c ) )
                adic = LOADFP( dintc( c ) )
                adix = LOADFP( dintx( c ) )
                adiv = LOADFP( ediv( c,b ) )

                b5 = FPADD( af1n, af4n )
                b5 = FPNMSUB( b5, two, af3n )
                b1 = FPNMSUB( b5, two, af2n )
                b2 = FPMADD( b1,three,af1n )
                b3 = FPNMSUB( af1t, two, af2t )
                b4 = FPADD( b3, af2n )
                b6 = FPNMSUB( b5, two, af4t )
                b4 = FPNMSUB( b4, two, af3t )
                b4 = FPNMSUB( b4, two, b6 )
                b4 = FXPMUL( b4, dintx(b) )
                b4 = FXCPMADD( b4, b1, dintc(b) )
                b4 = FPMUL( b4, at1v )
                b3 = FPMUL( b3, adic )
                b3 = FPMADD( b3, b2, adix )
                b3 = FXCPMADD( b4, b3, t1v(b) )
                axmp = FPMADD( axmp, b3, adiv ) 
            enddo
        enddo

        emp5 = DBLE( axmp ) + IMAG( axmp )

        bxmp = ( 0d0,0d0 )
        do b = 1, nvir, 2
            do c = 1, nvir, 2
                
                af1t = LOADFP( f1t( c,b ) )
                af2t = LOADFP( f2t( c,b ) )
                af3t = LOADFP( f3t( c,b ) )
                af4t = LOADFP( f4t( c,b ) )
                af1n = LOADFP( f1n( c,b ) )
                af2n = LOADFP( f2n( c,b ) )
                af3n = LOADFP( f3n( c,b ) )
                af4n = LOADFP( f4n( c,b ) )
                adiv = LOADFP( ediv( c,b ) )
                
                c1 = FPADD( af1t, af3n )
                c2 = FPADD( af1n, af2t )
                c2 = FPADD( c2, af4n )
                c9 = FPADD( af1t, af4t )
                c0 = FPADD( af2t, af3t )
                c3 = FPNMSUB( c9, two, c0 )
                c4 = FPADD( c1, af2n )
                c5 = FPMSUB( c0, two, c9 )
                c6 = FPMADD( af3n, two, af4t )
                c7 = FPMUL( af3n, af4t )
                c8 = FPMUL( c1, c3 )
                ca = FPMADD( c7, af1n, af1n )
                c8 = FPNMSUB( c8, c5, af1n )
                c8 = FPMADD( c8, three, ca )

                bf1t = LOADFP( f1t( b,c ) )
                bf2t = LOADFP( f2t( b,c ) )
                bf3t = LOADFP( f3t( b,c ) )
                bf4t = LOADFP( f4t( b,c ) )
                bf1n = LOADFP( f1n( b,c ) )
                bf2n = LOADFP( f2n( b,c ) )
                bf3n = LOADFP( f3n( b,c ) )

                b0 = FPADD( bf2t, bf3t )
                b9 = FPADD( bf1t, bf4t )
                b1 = FPADD( bf1t, bf3n )
                b3 = FPNMSUB( b9, two, b0 )
                b5 = FPMSUB( b0, two, b9 )

                cf1t = LOADFP( f1t( c,b+1 ) )
                cf2t = LOADFP( f2t( c,b+1 ) )
                cf3t = LOADFP( f3t( c,b+1 ) )
                cf4t = LOADFP( f4t( c,b+1 ) )
                cf1n = LOADFP( f1n( c,b+1 ) )
                cf2n = LOADFP( f2n( c,b+1 ) )
                cf3n = LOADFP( f3n( c,b+1 ) )
                cf4n = LOADFP( f4n( c,b+1 ) )
                bdiv = LOADFP( ediv( c,b+1 ) )

                d1 = FPADD( cf1t, cf3n )
                d2 = FPADD( cf1n, cf2t )
                d2 = FPADD( d2, cf4n )
                d9 = FPADD( cf1t, cf4t )
                d0 = FPADD( cf2t, cf3t )
                d3 = FPNMSUB( d9, two, d0 )
                d4 = FPADD( d1, cf2n )
                d5 = FPMSUB( d0, two, d9 )
                d6 = FPMADD( cf3n, two, cf4t )
                !d7 = FPMUL( cf3n, cf4t )
                !d8 = FPMUL( d1, d3 )
                !da = FPMADD( d7, cf1n, cf1n )
                !d8 = FPNMSUB( d8, d5, cf1n )
                !d8 = FPMADD( d8, three, da )
                d5 = FPNMSUB( d5, three, cf1n )
                d7 = FPMUL( d1, d3 )
                d8 = FPNMSUB( d7, d5, cf1n )

                ef1t = LOADFP( f1t( b,c+1 ) )
                ef2t = LOADFP( f2t( b,c+1 ) )
                ef3t = LOADFP( f3t( b,c+1 ) )
                ef4t = LOADFP( f4t( b,c+1 ) )
                ef1n = LOADFP( f1n( b,c+1 ) )
                ef2n = LOADFP( f2n( b,c+1 ) )
                ef3n = LOADFP( f3n( b,c+1 ) )

                e0 = FPADD( ef2t, ef3t )
                e9 = FPADD( ef1t, ef4t )
                e1 = FPADD( ef1t, ef3n )
                e3 = FPNMSUB( e9, two, e0 )
                e5 = FPMSUB( e0, two, e9 )

                f0 = IMAG( b3 )
                b3 = ( DBLE(b3),DBLE(e3) )
                e3 = ( f0,IMAG(e3) )
                f0 = IMAG( b5 )
                b5 = ( DBLE(b5),DBLE(e5) )
                e5 = ( f0,IMAG(e5) )
                f0 = IMAG( bf1n )
                bf1n = ( DBLE(bf1n),DBLE(ef1n) )
                ef1n = ( f0,IMAG(ef1n) )
                f0 = IMAG( bf2n )
                bf2n = ( DBLE(bf2n),DBLE(ef2n) )
                ef2n = ( f0,IMAG(ef2n) )

                c8 = FPMADD( c8, c2, b3 )
                c6 = FPMUL( c6, bf1n )
                c8 = FPNMSUB( c8, c4, b5 )
                c6 = FPMADD( c6, af2t, bf2n )
                c8 = FPMADD( c8, three, c6 )
                bxmp = FPMADD( bxmp, c8, adiv )

                d8 = FPMADD( d8, d2, e3 )
                d6 = FPMUL( d6, ef1n )
                d8 = FPNMSUB( d8, d4, e5 )
                d6 = FPMADD( d6, cf2t, ef2n )
                d6 = FPMADD( d6, cf3n, cf4t )
                d8 = FPMADD( d8, three, d6 )
                bxmp = FPMADD( bxmp, d8, bdiv )
            enddo
        enddo

        emp4 = emp4 + DBLE( bxmp ) + IMAG( bxmp )

        return
        end

c $Id$
