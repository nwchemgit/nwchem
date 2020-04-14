import math
import array

'''
A simple class and helper functions to read and compute with a
Guassian cube file as produced by NWChem with the DPLOT module.
It can read the file (or tabulate any function) and then compute the
value at an arbitrary interior point using tri-linear interpolation.

NWChem seems to produce the cube and atomic coordinates in atomic
units.

Written for Python 2.*.
'''

class Cube:
    '''
    3D uniform grid with bounding box that supports linear interpolation.

    Data members of the class are Nx, Ny, Nz (the number of points in
    each dimension), r0 and r1 the lower left and top right corners of
    the volume, dx, dy, dz the increments in each dimension.

    In the interpolation routine the volume axes are assumed aligned
    with the corresponding Cartesian axes.
    '''
    def __init__(self,N,r0,r1,f=None):
        '''
        Makes a Cube of dimension N=(Nx,Ny,Nz) with bounding corners r0,r1
        optionally initialized with function f(x,y,z).  If f is None, it
        is initialized to zero.
        '''
        
        self.d = array.array('d',[0.0]*N[0]*N[1]*N[2])
        self.Nx, self.Ny, self.Nz = N[0], N[1], N[2]
        self.r0, self.r1 = r0, r1
        self.dx, self.dy, self.dz = (r1[0]-r0[0])/(N[0]-1), (r1[1]-r0[1])/(N[1]-1), (r1[2]-r0[2])/(N[2]-1)

        # Used to fuzz numerical comparisons
        self.eps = 1e-15*max(self.dx*self.Nx,self.dy*self.Ny,self.dz*self.Nz)

        #print 'Cube:', self.r0, self.r1
        #print '   N:', self.Nx, self.Ny, self.Nz
        #print '   D:', self.dx, self.dy, self.dz
        #print ' eps:', self.eps
        if f:
            for iz in range(self.Nz):
                for iy in range(self.Ny):
                    for ix in range(self.Nx):
                        self[(ix,iy,iz)] = f(*self.get_coords(ix,iy,iz))

    def __getitem__(self,ind):
        ''' Gets/returns value indexed by index '''
        ix,iy,iz = ind
        return self.d[ix + self.Nx*(iy + self.Ny*iz)]

    def __setitem__(self,ind,value):
        ''' Assigns/sets value indexed by index '''
        ix,iy,iz = ind
        self.d[ix + self.Nx*(iy + self.Ny*iz)] = value

    def get_coords(self,ix,iy,iz):
        ''' Returns coordinates of index '''
        return (self.r0[0]+ix*self.dx,self.r0[1]+iy*self.dy,self.r0[2]+iz*self.dz)

    def containing_box(self,x,y,z):
        ''' Returns lower left index of box containing coords '''
        fac = 0.999999999999999
        ix, iy, iz = int((x-self.r0[0])*fac/self.dx), int((y-self.r0[1])*fac/self.dy), int((z-self.r0[2])*fac/self.dz)
        if ix<0 or ix>=(self.Nx-1) or iy<0 or iy>=(self.Ny-1) or iz<0 or iz>=(self.Nz-1):
            print("Trying to find containing box for point out of bounds")
            print("point", (x,y,z))
            print("lower", self.r0)
            print("upper", self.r1)
            raise IndexError
        return ix, iy, iz

    def linear_1d(self,x,xlo,xhi,flo,fhi):
        ''' Linear interpolation in 1D ... [xlo---x---xhi] '''

        if (xlo-x)>self.eps or (x-xhi)>self.eps:
            print("linear_1d: extrapolating!  (x, xlo, xhi):", x, xlo, xhi)
            raise IndexError
        return flo + (fhi-flo)*(x-xlo)/(xhi-xlo)
    
    def interp(self,x,y,z):
        ''' Linear interpolation within the 3D volume '''
        ix, iy, iz = self.containing_box(x,y,z)
        xlo,ylo,zlo = self.get_coords(ix,iy,iz)

        v000 = self[(ix  , iy  , iz  )]
        v100 = self[(ix+1, iy  , iz  )]
        v010 = self[(ix  , iy+1, iz  )]
        v001 = self[(ix  , iy  , iz+1)]
        v110 = self[(ix+1, iy+1, iz  )]
        v101 = self[(ix+1, iy  , iz+1)]
        v011 = self[(ix  , iy+1, iz+1)]
        v111 = self[(ix+1, iy+1, iz+1)]

        # interp in z
        v00 = self.linear_1d(z,zlo,zlo+self.dz,v000,v001)
        v10 = self.linear_1d(z,zlo,zlo+self.dz,v100,v101)
        v01 = self.linear_1d(z,zlo,zlo+self.dz,v010,v011)
        v11 = self.linear_1d(z,zlo,zlo+self.dz,v110,v111)

        # interp in y
        v0 = self.linear_1d(y,ylo,ylo+self.dy,v00,v01)
        v1 = self.linear_1d(y,ylo,ylo+self.dy,v10,v11)

        # interp in x
        return self.linear_1d(x,xlo,xlo+self.dx,v0,v1)
        
def read_i_f_f_f(f):
    ''' Read line containing integer and three floats '''
    line = f.readline().lstrip().rstrip().split()
    return int(line[0]), float(line[1]), float(line[2]), float(line[3])

def read_atom(f):
    ''' Read line from Guassian cube file containing atomic info '''
    line = f.readline().lstrip().rstrip().split()
    return int(line[0]), (float(line[2]), float(line[3]), float(line[4]))

def load_gaussian(filename):
    '''
    Returns tuple (cube,atoms) loaded from Gaussian cube file.

    An atom in the list is the tuple (Z,(X,Y,Z)).

    The volume axes are assumed aligned with the corresponding
    Cartesian axes.
    '''
    f = open(filename,'r')
    f.readline() # discard two comment lines
    f.readline()

    natoms, xlo, ylo, zlo = read_i_f_f_f(f)
    Nx, dx, junk, junk = read_i_f_f_f(f)
    Ny, junk, dy, junk = read_i_f_f_f(f)
    Nz, junk, junk, dz = read_i_f_f_f(f)
    
    # Load geometry
    atoms = []
    for i in range(natoms): 
        atoms.append(read_atom(f))

    xhi = xlo + (Nx-1)*dx
    yhi = ylo + (Ny-1)*dy
    zhi = zlo + (Nz-1)*dz
    
    c = Cube((Nx,Ny,Nz),(xlo,ylo,zlo),(xhi,yhi,zhi))
    
    maxval = 0.0
    for ix in range(Nx):
        for iy in range(Ny):
            line = f.readline().lstrip().rstrip().split()
            n = 0
            for iz in range(Nz):
                try:
                    value = float(line[n])
                    maxval = max(maxval,abs(value))
                    c[(ix,iy,iz)] = value
                except:
                    print(line)
                    print(n)
                    raise IndexError
                n = n + 1
                if n == 6:
                    line = f.readline().lstrip().rstrip().split()
                    n = 0

    #print "maxval",maxval
    return c, atoms

if __name__ == "__main__":

    def testfun(x,y,z):
        return math.sin(x+y*0.7+z*0.5)
        #return math.exp(-(x**2 + y**2 + z**2))
        #return 0.1*x + 0.3*y + 0.7*z
    
    def test():
        N = (50,60,70)
        r0 = (-0.5,-0.5,-0.5)
        r1 = (0.0,0.0,0.0)

        c = Cube(N,r0,r1,testfun)

        # print the interpolated function and error out along a line
        dx,dy,dz = (r1[0]-r0[0])/10, (r1[1]-r0[1])/10, (r1[2]-r0[2])/10,
        print("   x      y      z      exact     interp     error")
        print("------ ------ ------ ---------- ---------- ----------")
        for i in range(11):
            x,y,z = r0[0]+i*dx, r0[1]+i*dy, r0[2]+i*dz
            numeric = c.interp(x,y,z)
            exact = testfun(x,y,z)
            print("%6.2f %6.2f %6.2f %10.6f %10.6f %10.2e" % (x,y,z,exact,numeric,exact-numeric))
    
    test()
    
    #c, atoms = load_gaussian("dens.cube")
    #print atoms
    #print 'Cube:', c.r0, c.r1
    #print '   N:', c.Nx, c.Ny, c.Nz
    #print '   D:', c.dx, c.dy, c.dz
    #print ' eps:', c.eps
    
    
