##
## nw_spectrum
##
## Kenneth Lopata
## Last modified: 2013-03-08
##
## Python script for parsing NWChem output for TDDFT/vspec excitation
## energies, and optionally Lorentzian broadenening the spectra.  For
## online help run:
##
## nw_spectrum --help
##
##

import sys
import textwrap
from optparse import OptionParser

ver = "2.1"
pname = "nw_spectrum"


def check_version ():
    """Check version and ensure new print and string formatting is
    allowed.  Raises and exception if not satisfied."""

    if sys.version_info < (2, 6):
        raise Exception("This script requires python >= 2.6")

    try:
        newstring = "oldstring {v}".format(v=3.14)
    except:
        raise Exception("This script requires string.format()")


def ev2au(e_ev):
    return (1.0 / 27.2114) * e_ev


def au2ev (e_au):
    return 27.2114 * e_au


def ev2nm(e_ev):
    return 27.2114 * 2.0 * 2.99792 * 2.41888 * 3.14159265359 / e_ev


def determine_data_type ():
    """Parses stdin to see what data type, then rewinds stdin.
    Returns 'vspec', 'tddft', or raises and Exception if neither
    found.  It choses based on the first tag found in the input, so
    files with multiple data will only find the 1st.  To extract the
    2nd, manually specify the data format via the command line args."""
    
    tag_tddft = "NWChem TDDFT Module"
    tag_vspec = "DFT Virtual Spectrum"

    lines = sys.stdin.readlines()
    
    for line in lines:
        if tag_tddft in line:
            sys.stdin.seek(0)
            return "tddft"
        elif tag_vspec in line:
            sys.stdin.seek(0)
            return "vspec"

    raise Exception ("Failed to determine data format, please specify manually.")


def parse_input_vspec (opts):
    """Parses input from vspec and returns excitation energies in the
    form [energy, f], in eV and atomic units units, respectively."""

    lines = sys.stdin.readlines ()
    
    inside_data = False
    roots = []
    for l in lines:
        
        if "<START>" in l:
            try:
                ls = l.split()
                tag = ls[0]
                nexcite = int (ls[1])
            except:
                raise Exception ("Failed to parse <START> tag and number: {0}".format(l))

            iexcite = 0
            inside_data = True
            continue

        if inside_data:
            if "<END>" in l:
                inside_data = False
                continue
#                break
            
            try:
                line_split = l.strip().split()
                n = int (line_split[0])
                occ = int (line_split[1])          #not used
                virtual = int (line_split[2])      #not used
                energy_ev = float (line_split[3])
                osc = float (line_split[7])
            except:
                raise Exception ("Failed to parse data line: {0}".format(l))
            
            iexcite = iexcite + 1
            
            if n != iexcite:
                raise Exception ("Expected excitation number {0}, found {1}".format(iexcite, n))
            
            
            if energy_ev < 0.0:
                print ("{0} Warning: Ignored negative vpsec excitation: {1} eV, {2}".format(opts.cchar, energy_ev, osc))
                
                if opts.verbose:
                    sys.stderr.write ("Warning: Ignored negative vpsec excitation: {0} eV, {1}\n".format(energy_ev, osc))
            else:
                roots.append ([energy_ev, osc])


#    if not inside_data:
#        raise Exception ("Failed to find <START> tag")

    if iexcite != nexcite:
        print ("{0} Warning: Expected {1} excitations, found {2}".format(opts.cchar, nexcite, iexcite))
        
        if opts.verbose:
            sys.stderr.write ("Warning: Expected {0} excitations, found {1}\n".format(nexcite,iexcite))


    if opts.verbose:
        sys.stderr.write ("{0}: Found {1} vspec excitations\n".format(pname, len(roots)))
            
            
    return roots


def parse_input_evals (opts):
    """Parses input for eigenvalues and return as a list"""

    start_tag = "DFT Final Molecular Orbital Analysis"
    end_tag = "Task  times  cpu"

    inside = False
    
    lines = sys.stdin.readlines ()
    iline = -1
    evals = []

    while True:
        iline = iline + 1
        try:
            line = lines[iline]
        except:
            break

        if start_tag in line:
            inside = True

        if end_tag in line:
            inside = False

        if inside and "Vector" in line and "Occ" in line:
            line_strip = line.strip()

            try:
                tagloc = line_strip.rfind("E=")
                evalue_str = line_strip[tagloc+2:tagloc+15].replace("D", "E")  # after E=  ; replace D with E
                evalue = float(evalue_str)
            except:
                raise Exception ("Failed to parse eigenvalue: {0}".format(line_strip))

            eval_ev = au2ev (evalue)
            evals.append(eval_ev)  #store eigenvalue in eV


    if opts.verbose:
        sys.stderr.write ("{0}: Found {1} eigenvalues\n".format(pname, len(evals)))

    return evals


def bin_evals (opts, evals):
    """Take eigenvalues and bins them, and return [energy, N], where N
    is the number of eigenvalues in the energy bin centered around
    energy"""

    ## they should be sorted but let's make sure
    evals_sort = sorted(evals)

    emin = evals_sort[0]
    emax = evals_sort[-1]

    de = (emax - emin) / opts.nbin
    
    #=== XXX HARDCODE ===
#    de = 0.01
#    opts.nbin = int((emax - emin)/de) + 1
    #=== 

    dos_raw = []
    for ie in range(opts.nbin+1):
        ecenter = emin + ie*de
        eleft = ecenter - 0.5*de
        eright = ecenter + 0.5*de

        count = 0
        for val in evals_sort:
            if val >= eleft and val <= eright:
                count = count + 1

        dos_raw.append ([ecenter, count])


    ## check that total sum is number of eigenvalues
    ntot = 0
    for d in dos_raw:
        ntot = ntot + d[1]

    if ntot != len (evals):
        raise Exception ("Inconsistent integrated DOS and number of eigenvalues: {0} vs {1}".format(ntot, len(evals)))

    return dos_raw


def parse_input_tddft (opts):
    """Parses input for singlet TDDFT excitation energies.
    Returns excitation energies in the form [energy, f], in eV and
    atomic units units, respectively."""

    start_tag = "Convergence criterion met"
    end_tag = "Excited state energy"
    max_osc_search = 10  #max number of lines after root energy to look for oscillator strength


    inside = False  #true when we are inside output block
    
    lines = sys.stdin.readlines()

    iline = -1
    roots = []
    
    while True:
        
        iline = iline + 1
        try:
            line = lines[iline]
        except:
            break

        if start_tag in line:
            inside = True

        if end_tag in line:
            inside = False

        if inside and "Root" in line and "eV" in line:

            line_strip = line.strip()

            # ##
            # ## OLD WAY
            # ## 
            # ## Note, I would ideally .split() the line directly, but you
            # ## can have cases like Root356 (for three digit root numbers)
            # ## so I have to do it this ugly way...
            # ##
            # try:
            #     line_start = line_strip[0:4]          # should contain RootXXX
            #     line_data = line_strip[7:].split()   
            #     line_n = line_strip[4:7]              #contains integer root number
            #     line_e = line_data[-2]                # should contain excitation energy
            #     line_ev_tag = line_data[-1]           # should contain "eV"
            # except:
            #     raise Exception ("Failed to parse data line for root: {0}".format(line_strip))

            ##
            ## NEW WAY after Niri changed formatting, e.g.:
            ## Root  48 singlet a             57.454053695 a.u.             1563.4050 eV
            ##
            line_strip = line.strip()
            line_split = line_strip.split()
            try:
                line_start = line_split[0]     # contains "Root"
                line_n = line_split[1]         # contains root number (int)
                line_ev_tag = line_split[-1]   # contains "eV"
                line_e = line_split[-2]        # contains excitation energy in eV
            except:
                raise Exception ("Failed to parse data line for root: {0}".format(line_strip))

            if line_start == "Root" and line_ev_tag == "eV":
                try:
                    n = int(line_n)
                    energy_ev = float(line_e)
                except:
                    raise Exception ("Failed to convert root values: {0}".format(line_strip))
            else:
                raise Exception ("Unexpected format for root: {0}".format(line_strip))

            
            if line_start == "Root" and line_ev_tag == "eV":
                try:
                    n = int(line_n)
                    energy_ev = float(line_e)
                except:
                    raise Exception ("Failed to convert root values: {0}".format(line_strip))
            else:
                raise Exception ("Unexpected format for root: {0}".format(line_strip))


            ##
            ## Now look for oscillator strength, which will be a few
            ## lines down (though the exact position may vary it seems).
            ##
            ioscline = -1
            while True:
                ioscline = ioscline + 1
                if ioscline >= max_osc_search:
                    raise Exception ("Failed to find oscillator strength after looking {0} lines.".format(ioscline))
                
                oscline = lines[iline + ioscline].strip()
                
                if "Dipole Oscillator Strength" in oscline:
                    try:
                        osc_str = oscline.split()
                        osc = float (osc_str[3])
                    except:
                        raise Exception ("Failed to convert oscillator strength: {0}".format(oscline))
                    break
            
            
            ## do some final checks, then append to data
            if energy_ev < 0.0:
                raise Exception ("Invalid negative energy: {0}".format(energy_ev))

            if osc < 0.0:
                raise Exception ("Invalid negative oscillator strength: {0}".format(osc))

            roots.append([energy_ev, osc])


    nroots = len (roots)

    if nroots < 1:
        raise Exception ("Failed to find any TDDFT roots")
    else:
        if opts.header:
            print ("{0} Successfully parsed {1} TDDFT singlet excitations".format(opts.cchar,nroots))


    if opts.verbose:
        sys.stderr.write ("{0}: Found {1} TDDFT excitations\n".format(pname, nroots))

    return roots



def make_energy_list (opts, roots):
    """Computes the list of spectrum energies, and potentially adjusts
    peak widths"""

    epad = 20.0*opts.width
    
    emin = roots[0][0] - epad

#    if emin < opts.width:
#        emin = opts.width

    emax = roots[-1][0] + epad
    de = (emax - emin) / opts.npoints

    ## Use width of at least two grid points
    if opts.width < 2*de:
        opts.width = 2*de
        print ("{0} Warning: Forced broadening to be {1} eV".format(opts.cchar, opts.width))

        if opts.verbose:
            sys.stderr.write ("Warning: Forced broadening to be {0} eV\n".format(opts.width))
    
#    opts.width = max (opts.width, 2*de)

    eout = [ emin + ie*de for ie in range(opts.npoints) ]
    return eout


### OLD SLOWER WAY
# def lorentzian_broaden (opts, roots):
#     """Broadens raw roots into spectrum and returns the result."""

#     ## cutoff radius
#     cutoff = 15.0*opts.width

#     ##
#     ## multiply by 0.5 as FWHM was supplied by gamma in lorenzian is
#     ## actually HWHM.
#     ##
#     gamma = 0.5 * opts.width
#     gamma_sqrd = gamma*gamma


#     ##
#     ## L(w; w0, gamma) = gamma/pi * 1 / [(w-w0)^2 + gamma^2]
#     ##
#     prefac = gamma/3.14159265359*de


#     # spectrum = [ [emin + ie*de, 0] for ie in range(opts.npoints)]
    
#     # for point in spectrum:
#     #     stot = 0.0
#     #     for root in roots:
#     #         xx0 = point[0] - root[0]
        
#     #         if abs(xx0) <= cutoff:
#     #             stot += prefac * root[1] / ( xx0*xx0 + gamma_sqrd)   #Lorentzian
            
#     #         point[1] = stot


#     # npoints_made = len(spectrum)
#     # if npoints_made != opts.npoints:
#     #     raise Exception ("Spectrum should have {0} points, instead has {1}".format(opts.npoints, npoints_made))

#     # if opts.header:
#     #     print ("{0} Lorentzian broadened roots with width {1} eV".format(opts.cchar, opts.width))
#     #     print ("{0} Created spectrum with {1} points".format(opts.cchar, npoints_made))

#     return spectrum



## Faster way--use a generator
def gen_spectrum (opts, energies, roots):
    """Generator for making Lorentzian broadenend spectrum."""
    
    ## cutoff radius
#    cutoff = 15.0*opts.width
    cutoff = 20.0*opts.width
#    cutoff = 9999999.0*opts.width

    
    ##
    ## L(w; w0, gamma) = gamma/pi * 1 / [(w-w0)^2 + gamma^2]
    ##
    ##
    ## multiply by 0.5 as FWHM was supplied by gamma in
    ## lorenzian is actually HWHM.
    ##
    gamma = 0.5 * opts.width
    gamma_sqrd = gamma*gamma

    de = (energies[-1] - energies[0]) / (len(energies)-1)
    prefac = gamma/3.14159265359*de

    
    ##
    ## used for verbose output (ie see progress for very large runs)
    ##
    # checkpt_energies = []
    # ne = len(energies)
    # for ie in range(ne):
    #     if ie % ne == 0:
    #         checkpt_energies.append(energies[ie])


    for energy in energies:

        # if opts.verbose:
        #     if energy in checkpt_energies:
        #         sys.stderr.write ("XXX\n")
            
        
        stot = 0.0
        for root in roots:
            xx0 = energy - root[0]
        
            if abs(xx0) <= cutoff:
                stot += root[1] / ( xx0*xx0 + gamma_sqrd)   #Lorentzian
            
        yield [energy, stot*prefac]



        

def dump_header (opts):
    """Print header to stdout"""

    if opts.header:
        print ("{0} ================================".format(opts.cchar))
        print ("{0}  NWChem spectrum parser ver {1}".format(opts.cchar,ver))
        print ("{0} ================================".format(opts.cchar))
        print ("{0} ".format(opts.cchar))
        print ("{0} Parser runtime options: {1}".format(opts.cchar,opts))
        print ("{0}".format(opts.cchar))


def dump_data (opts, roots):
    """Dumps output to stdout.  This works for either lists of raw roots
    or broadened spectra."""

    if opts.verbose:
        sys.stderr.write ("{0}: Dumping data to stdout ... \n".format(pname))

    if opts.units == "ev":
        if opts.header:
            print ("{0}".format(opts.cchar))
            print ("{c}{s1}{delim}{s2}".format(c=opts.cchar, s1="  Energy [eV]  ", delim=opts.delim, s2="   Abs. [au]   "))
            print ("{0}----------------------------------".format(opts.cchar))
            
        for root in roots:
            print ("{energy:.10e}{delim}{osc:.10e}".format(energy=root[0], delim=opts.delim, osc=root[1]))


    elif opts.units == "au":
        if opts.header:
            print ("{0}".format(opts.cchar))
            print ("{c}{s1}{delim}{s2}".format(c=opts.cchar, s1="  Energy [au]  ", delim=opts.delim, s2="   Abs. [au]   "))
            print ("{0}----------------------------------".format(opts.cchar))
            
        for root in roots:
            eout = ev2au (root[0])
            print ("{energy:.10e}{delim}{osc:.10e}".format(energy=eout, delim=opts.delim, osc=root[1]))

    elif opts.units == "nm":
        if opts.header:
            print ("{0}".format(opts.cchar))
            print ("{c}{s1}{delim}{s2}".format(c=opts.cchar, s1=" Wavelen. [nm] ", delim=opts.delim, s2="   Abs. [au]   "))
            print ("{0}----------------------------------".format(opts.cchar))
            
#        roots.reverse ()  #reorder so we output in increasing wavelength
#XXX SHOULD REVERSE            
        for root in roots:
            eout = ev2nm (root[0])
            print ("{energy:.10e}{delim}{osc:.10e}".format(energy=eout, delim=opts.delim, osc=root[1]))

    else:
        raise Exception ("Invalid unit: {0}".format(opts.units))        
        

def preprocess_check_opts (opts):
    """Check options and replaces with sane values if needed, stores
    all opts as purely lowercase."""

    opts.units = opts.units.lower()
    opts.datafmt = opts.datafmt.lower()

    if opts.units != "nm" and opts.units != "ev" and opts.units != "au":
        raise Exception ("Invalid unit type: {0}".format(opts.units))

    if opts.datafmt != "tddft" and opts.datafmt != "vspec" and opts.datafmt != "auto" and opts.datafmt != "dos":
        raise Exception ("Invalid data format: {0}".format(opts.datafmt))
        
    if opts.npoints < 100:
        raise Exception ("Increase number of points to at least 100 (you asked for {0})".format(opts.npoints))

    if opts.width < 0.0:
        raise Exception ("Peak width must be positive (you supplied {0})".format(opts.width))


def main():
    
    ##
    ## Parse command line options.  Note we later make all lowercase,
    ## so from the user's point of view they are case-insensitive.
    ##
    usage = "%prog [options]\n\n"
    desc = "Reads NWChem output from stdin, parses for the linear response TDDFT or DFT vspec excitations, and prints the absorption spectrum to stdout.  It will optionally broaden peaks using a Lorentzian with FWHM of at least two energy/wavelength spacings.  By default, it will automatically determine data format (tddft or vspec) and generate a broadened spectrum in eV."
    desc_wrap = textwrap.wrap (desc,80)
    for s in desc_wrap:
        usage += s + "\n"

    example = 'Create absorption spectrum in nm named "spectrum.dat" from the NWChem output file "water.nwo" named spectrum.dat with peaks broadened by 0.3 eV and 5000 points in the spectrum.'
    example_wrap = textwrap.wrap (example,80)
    usage += "\nExample:\n\n\t"+"nw_spectrum -b0.3 -p5000 -wnm < water.nwo > spectrum.dat\n\n"
    for s in example_wrap:
        usage += s + "\n"


    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--format", type="string", dest="datafmt",
                      help="data file format: auto (default), tddft, vspec, dos", metavar="FMT")
    parser.add_option("-b", "--broad", type="float", dest="width", 
                      help="broaden peaks (FWHM) by WID eV (default 0.1 eV)", metavar="WID")
    parser.add_option("-n", "--nbin", type="int", dest="nbin",
                      help="number of eigenvalue bins for DOS calc (default 20)", metavar="NUM")
    parser.add_option("-p", "--points", type="int", dest="npoints",
                      help="create a spectrum with NUM points (default 2000)", metavar="NUM")
    parser.add_option("-w", "--units", type="string", dest="units",
                      help="units for frequency:  eV (default), au, nm", metavar="UNT")
    parser.add_option("-d", "--delim", type="string", dest="delim",
                      help="use STR as output separator (four spaces default)", metavar="STR")
    parser.add_option("-x", "--extract", action="store_false", dest="makespec",
                      help="extract unbroadened roots; do not make spectrum")
    parser.add_option("-C", "--clean", action="store_false", dest="header",
                      help="clean output; data only, no header or comments")
    parser.add_option("-c", "--comment", type="string", dest="cchar",
                      help="comment character for output ('#' default)", metavar="CHA")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="echo warnings and progress to stderr")


    parser.set_defaults(width = 0.1, npoints = 2000, units="eV", datafmt="auto", delim="    ", makespec=True, header=True, cchar="#", nbin=20, verbose=False)
    (opts, args) = parser.parse_args()

    preprocess_check_opts(opts)
    check_version ()
    
    if opts.header:
        dump_header (opts)

    if opts.datafmt == "auto":
        opts.datafmt = determine_data_type ()
        if opts.datafmt == "tddft":
            print ("{0} The input appears to contain TDDFT data.".format(opts.cchar))
        elif opts.datafmt == "vspec":
            print ("{0} The input appears to contain vspec data.".format(opts.cchar))
        else:
            raise Exception ("Invalid data format: {0}".format(opts.datafmt))
            

    ## parse raw data
    if opts.datafmt == "tddft":
        roots = parse_input_tddft (opts)
    elif opts.datafmt == "vspec":
        roots = parse_input_vspec (opts)
    elif opts.datafmt == "dos":
        evals = parse_input_evals (opts)
        roots = bin_evals (opts, evals)
    else:
        raise Exception ("Invalid data format supplied: {0}".format(opts.datafmt))

    
    ## make broadened spectrum if desired
    if opts.makespec:
        energies = make_energy_list (opts, roots)
        spectrum = gen_spectrum (opts, energies, roots)

        if opts.verbose:
            sys.stderr.write ("{0}: Initialized spectrum [{1: .3f} : {2: .3f}] ({3} points)\n".format(pname, energies[0], energies[-1], len(energies)))

        if opts.header:
            print ("{0} Roots were Lorentzian broadened with width {1} eV ".format(opts.cchar, opts.width))
            print ("{0} Spectrum generated with {1} points.".format(opts.cchar, opts.npoints))
            
    else:
        spectrum = roots
        if opts.header:
            print ("{0} Roots not broadened".format(opts.cchar))
            print ("{0} No spectrum generated: output is list of raw excitations".format(opts.cchar))


    ## finally, dump data to stdout
    dump_data (opts, spectrum)
        
    if opts.verbose:
        sys.stderr.write ("{0}: Done\n".format(pname))


if __name__ == "__main__":
    main()


