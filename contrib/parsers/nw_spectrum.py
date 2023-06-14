#!/usr/bin/env python
##
## nw_spectrum
##
## Kenneth Lopata
## Last modified: 2015-02-21
##
## Python script for parsing NWChem output for TDDFT/vspec excitation
## energies, and optionally Lorentzian broadenening the spectra.  For
## online help run:
##
## nw_spectrum.py --help
##
##
## Examples, using data from NWChem QA reference output:
## tddft: nw_spectrum.py -v -i QA/tests/tddft_h2o/tddft_h2o.out -o tddft.dat
## vspec: TODO, no test cases exist!
## dos: nw_spectrum.py -v -f dos -i QA/tests/tddft_h2o/tddft_h2o.out -o dos.dat
##
## To plot an example file, in gnuplot, do e.g.
## plot "tddft.dat" with lines

import sys
import textwrap
from optparse import OptionParser

ver = "2.2"
pname = "nw_spectrum"


class SpectrumProcessor(object):
    def __init__(self, opts, instream, outstream):
        self.instream = instream
        self.outstream = outstream
        self.opts = opts
        self.preprocess_check_opts()
        self.check_version()

    def write(self, line):
        """Write a line of output to the output stream"""

        self.outstream.write(line + "\n")
        
    def check_version(self):
        """Check version and ensure new print and string formatting is
        allowed.  Raises and exception if not satisfied."""

        if sys.version_info < (2, 6):
            raise Exception("This script requires python >= 2.6")

        try:
            newstring = "oldstring {v}".format(v=3.14)
        except:
            raise Exception("This script requires string.format()")


    def ev2au(self, e_ev):
        return (1.0 / 27.2114) * e_ev

    def au2ev(self, e_au):
        return 27.2114 * e_au

    def ev2nm(self, e_ev):
        return 27.2114 * 2.0 * 2.99792 * 2.41888 * 3.14159265359 / e_ev

    def determine_data_type(self):
        """Parses input to see what data type, then rewinds input.
        Returns 'vspec' or 'tddft', or raises an exception if neither
        found.  It chosses based on the first tag found in the input, so
        files with multiple data will only find the 1st.  To extract the
        2nd, manually specify the data format via the command line args."""

        tag_tddft = "NWChem TDDFT Module"
        tag_vspec = "DFT Virtual Spectrum"

        lines = self.instream.readlines()

        for line in lines:
            if tag_tddft in line:
                self.instream.seek(0)
                return "tddft"
            elif tag_vspec in line:
                self.instream.seek(0)
                return "vspec"

        raise Exception("Failed to determine data format, please specify manually.")


    def parse_input_vspec(self):
        """Parses input from vspec and returns excitation energies in the
        form [energy, f], in eV and atomic units units, respectively."""

        lines = self.instream.readlines()

        inside_data = False
        roots = []
        for l in lines:
            if "<START>" in l:
                try:
                    ls = l.split()
                    tag = ls[0]
                    nexcite = int (ls[1])
                except:
                    raise Exception("Failed to parse <START> tag and number: {0}".format(l))

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
                    raise Exception("Failed to parse data line: {0}".format(l))

                iexcite = iexcite + 1

                if n != iexcite:
                    raise Exception("Expected excitation number {0}, found {1}".format(iexcite, n))

                if energy_ev < 0.0:
                    self.write("{0} Warning: Ignored negative vpsec excitation: {1} eV, {2}".format(opts.cchar, energy_ev, osc))
                    if self.opts.verbose:
                        sys.stderr.write ("Warning: Ignored negative vpsec excitation: {0} eV, {1}\n".format(energy_ev, osc))
                else:
                    roots.append([energy_ev, osc])


    #    if not inside_data:
    #        raise Exception("Failed to find <START> tag")

        if iexcite != nexcite:
            self.write("{0} Warning: Expected {1} excitations, found {2}".format(self.opts.cchar, nexcite, iexcite))

            if self.opts.verbose:
                sys.stderr.write("Warning: Expected {0} excitations, found {1}\n".format(nexcite, iexcite))


        if self.opts.verbose:
            sys.stderr.write("{0}: Found {1} vspec excitations\n".format(pname, len(roots)))

        return roots


    def parse_input_evals(self):
        """Parses input for eigenvalues and return as a list"""

        start_tag = "DFT Final Molecular Orbital Analysis"
        end_tag = "Task  times  cpu"

        inside = False

        lines = self.instream.readlines()
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
                    evalue_str = line_strip[tagloc + 2 : tagloc + 15].replace("D", "E")  # after E=  ; replace D with E
                    evalue = float(evalue_str)
                except:
                    raise Exception("Failed to parse eigenvalue: {0}".format(line_strip))

                eval_ev = self.au2ev(evalue)
                evals.append(eval_ev)  #store eigenvalue in eV

        if self.opts.verbose:
            sys.stderr.write("{0}: Found {1} eigenvalues\n".format(pname, len(evals)))

        return evals


    def bin_evals(self, evals):
        """Take eigenvalues and bins them, and return [energy, N], where N
        is the number of eigenvalues in the energy bin centered around
        energy"""

        ## they should be sorted but let's make sure
        evals_sort = sorted(evals)

        emin = evals_sort[0]
        emax = evals_sort[-1]

        de = (emax - emin) / self.opts.nbin

        #=== XXX HARDCODE ===
    #    de = 0.01
    #    self.opts.nbin = int((emax - emin)/de) + 1
        #=== 

        dos_raw = []
        for ie in range(self.opts.nbin + 1):
            ecenter = emin + ie * de
            eleft = ecenter - 0.5 * de
            eright = ecenter + 0.5 * de

            count = 0
            for val in evals_sort:
                if val >= eleft and val <= eright:
                    count += 1

            dos_raw.append([ecenter, count])


        ## check that total sum is number of eigenvalues
        ntot = 0
        for d in dos_raw:
            ntot = ntot + d[1]

        if ntot != len(evals):
            raise Exception("Inconsistent integrated DOS and number of eigenvalues: {0} vs {1}".format(ntot, len(evals)))

        return dos_raw


    def parse_input_tddft(self):
        """Parses input for singlet TDDFT excitation energies.
        Returns excitation energies in the form [energy, f], in eV and
        atomic units units, respectively."""

        start_tag = "Convergence criterion met"
        end_tag = "Excited state energy"
        singlet_tag = "singlet excited"
        triplet_tag = "triplet excited"
        state = "singlet"

        max_osc_search = 10  #max number of lines after root energy to look for oscillator strength


        inside = False  #true when we are inside output block

        lines = self.instream.readlines()

        iline = -1
        roots = []

        while True:
            iline = iline + 1
            try:
                line = lines[iline]
            except:
                break

            #Only extract singlet state data. By default TDDFT calculations
            #produce both singlet and triplet data, and the (possibly
            #meaningless) triplet state data can otherwise interfere with
            #proper parsing of singlets.
            if start_tag in line and state == "singlet":
                inside = True

            if end_tag in line:
                inside = False

            if singlet_tag in line:
                state = "singlet"

            if triplet_tag in line:
                state = "triplet"

            if inside and "Root" in line and "eV" in line:
                line_strip = line.strip()
                line_split = line_strip.split()
                try:
                    line_start = line_split[0]     # contains "Root"
                    line_n = line_split[1]         # contains root number (int)
                    line_ev_tag = line_split[-1]   # contains "eV"
                    line_e = line_split[-2]        # contains excitation energy in eV
                except:
                    raise Exception("Failed to parse data line for root: {0}".format(line_strip))

                if line_start == "Root" and line_ev_tag == "eV":
                    try:
                        n = int(line_n)
                        energy_ev = float(line_e)
                    except:
                        raise Exception("Failed to convert root values: {0}".format(line_strip))
                else:
                    raise Exception("Unexpected format for root: {0}".format(line_strip))

                if line_start == "Root" and line_ev_tag == "eV":
                    try:
                        n = int(line_n)
                        energy_ev = float(line_e)
                    except:
                        raise Exception("Failed to convert root values: {0}".format(line_strip))
                else:
                    raise Exception("Unexpected format for root: {0}".format(line_strip))
                ##
                ## Now look for oscillator strength, which will be a few
                ## lines down (though the exact position may vary it seems).
                ##
                ioscline = -1
                while True:
                    ioscline = ioscline + 1
                    if ioscline >= max_osc_search:
                        raise Exception("Failed to find oscillator strength after looking {0} lines.".format(ioscline))

                    oscline = lines[iline + ioscline].strip()

                    if "Dipole Oscillator Strength" in oscline:
                        try:
                            osc_str = oscline.split()
                            osc = float(osc_str[3])
                        except:
                            raise Exception("Failed to convert oscillator strength: {0}".format(oscline))
                        break

                ## do some final checks, then append to data
                if energy_ev < 0.0:
                    raise Exception("Invalid negative energy: {0}".format(energy_ev))

                if osc < 0.0:
                    raise Exception("Invalid negative oscillator strength: {0}".format(osc))

                roots.append([energy_ev, osc])


        nroots = len(roots)

        if nroots < 1:
            raise Exception("Failed to find any TDDFT roots")
        else:
            if self.opts.header:
                self.write("{0} Successfully parsed {1} TDDFT singlet excitations".format(self.opts.cchar, nroots))


        if self.opts.verbose:
            sys.stderr.write("{0}: Found {1} TDDFT excitations\n".format(pname, nroots))

        return roots


    def make_energy_list(self, roots):
        """Computes the list of spectrum energies, and potentially adjusts
        peak widths"""

        epad = 20.0 * self.opts.width

        emin = roots[0][0] - epad

        #if emin < self.opts.width:
        #    emin = self.opts.width

        emax = roots[-1][0] + epad
        de = (emax - emin) / self.opts.npoints

        ## Use width of at least two grid points
        if self.opts.width < 2 * de:
            self.opts.width = 2 * de
            if self.opts.verbose:
                sys.stderr.write("Warning: Forced broadening to be {0} eV\n".format(self.opts.width))

        eout = [emin + ie * de for ie in range(self.opts.npoints)]
        return eout

    def gen_spectrum(self, energies, roots):
        """Generator for making Lorentzian broadenend spectrum."""

        ## cutoff radius
        #cutoff = 15.0 * self.opts.width
        cutoff = 20.0 * self.opts.width
        #cutoff = 9999999.0 * self.opts.width


        ##
        ## L(w; w0, gamma) = gamma/pi * 1 / [(w-w0)^2 + gamma^2]
        ##
        ##
        ## multiply by 0.5 as FWHM was supplied by gamma in
        ## lorenzian is actually HWHM.
        ##
        gamma = 0.5 * self.opts.width
        gamma_sqrd = gamma * gamma

        de = (energies[-1] - energies[0]) / (len(energies) - 1)
        prefac = gamma / 3.14159265359 * de


        ##
        ## used for verbose output (ie see progress for very large runs)
        ##
        # checkpt_energies = []
        # ne = len(energies)
        # for ie in range(ne):
        #     if ie % ne == 0:
        #         checkpt_energies.append(energies[ie])


        for energy in energies:
            # if self.opts.verbose:
            #     if energy in checkpt_energies:
            #         sys.stderr.write("XXX\n")


            stot = 0.0
            for root in roots:
                xx0 = energy - root[0]

                if abs(xx0) <= cutoff:
                    stot += root[1] / ( xx0 * xx0 + gamma_sqrd)   #Lorentzian

            yield [energy, stot * prefac]

    def dump_header(self):
        """Add header to output"""

        if self.opts.header:
            self.write("{0} ================================".format(self.opts.cchar))
            self.write("{0}  NWChem spectrum parser ver {1}".format(self.opts.cchar, ver))
            self.write("{0} ================================".format(self.opts.cchar))
            self.write("{0} ".format(self.opts.cchar))
            self.write("{0} Parser runtime options: {1}".format(self.opts.cchar, self.opts))
            self.write("{0}".format(self.opts.cchar))


    def dump_data(self, roots):
        """Dumps output.  This works for either lists of raw roots or broadened
        spectra."""

        if self.opts.verbose:
            sys.stderr.write("{0}: Dumping data ... \n".format(pname))

        if self.opts.units == "ev":
            if self.opts.header:
                self.write("{0}".format(self.opts.cchar))
                self.write("{c}{s1}{delim}{s2}".format(c=self.opts.cchar, s1="  Energy [eV]  ", delim=self.opts.delim, s2="   Abs. [au]   "))
                self.write("{0}----------------------------------".format(self.opts.cchar))

            for root in roots:
                self.write("{energy:.10e}{delim}{osc:.10e}".format(energy=root[0], delim=self.opts.delim, osc=root[1]))


        elif self.opts.units == "au":
            if self.opts.header:
                self.write("{0}".format(self.opts.cchar))
                self.write("{c}{s1}{delim}{s2}".format(c=self.opts.cchar, s1="  Energy [au]  ", delim=self.opts.delim, s2="   Abs. [au]   "))
                self.write("{0}----------------------------------".format(self.opts.cchar))

            for root in roots:
                eout = self.ev2au(root[0])
                self.write("{energy:.10e}{delim}{osc:.10e}".format(energy=eout, delim=self.opts.delim, osc=root[1]))

        elif self.opts.units == "nm":
            if self.opts.header:
                self.write("{0}".format(self.opts.cchar))
                self.write("{c}{s1}{delim}{s2}".format(c=self.opts.cchar, s1=" Wavelen. [nm] ", delim=self.opts.delim, s2="   Abs. [au]   "))
                self.write("{0}----------------------------------".format(self.opts.cchar))

    #        roots.reverse ()  #reorder so we output in increasing wavelength
    #XXX SHOULD REVERSE
            for root in roots:
                eout = self.ev2nm(root[0])
                self.write("{energy:.10e}{delim}{osc:.10e}".format(energy=eout, delim=self.opts.delim, osc=root[1]))

        else:
            raise Exception("Invalid unit: {0}".format(self.opts.units))

    def preprocess_check_opts(self):
        """Check options and replaces with sane values if needed, stores
        all opts as purely lowercase."""

        self.opts.units = self.opts.units.lower()
        self.opts.datafmt = self.opts.datafmt.lower()

        if self.opts.units not in ("nm", "ev", "au"):
            raise Exception("Invalid unit type: {0}".format(self.opts.units))

        if self.opts.datafmt not in ("tddft", "vspec", "auto", "dos"):
            raise Exception("Invalid data format: {0}".format(self.opts.datafmt))

        if self.opts.npoints < 100:
            raise Exception("Increase number of points to at least 100 (you asked for {0})".format(self.opts.npoints))

        if self.opts.width < 0.0:
            raise Exception("Peak width must be positive (you supplied {0})".format(self.opts.width))

    def process_all(self):
        if self.opts.header:
            self.dump_header()

        if self.opts.datafmt == "auto":
            self.opts.datafmt = self.determine_data_type()
            if self.opts.datafmt == "tddft":
                self.write("{0} The input appears to contain TDDFT data.".format(self.opts.cchar))
            elif self.opts.datafmt == "vspec":
                self.write("{0} The input appears to contain vspec data.".format(self.opts.cchar))
            else:
                raise Exception("Invalid data format: {0}".format(self.opts.datafmt))

        ## parse raw data
        if self.opts.datafmt == "tddft":
            roots = self.parse_input_tddft()
        elif self.opts.datafmt == "vspec":
            roots = self.parse_input_vspec()
        elif self.opts.datafmt == "dos":
            evals = self.parse_input_evals()
            roots = self.bin_evals(evals)
        else:
            raise Exception("Invalid data format supplied: {0}".format(self.opts.datafmt))


        ## make broadened spectrum if desired
        if self.opts.makespec:
            energies = self.make_energy_list(roots)
            spectrum = self.gen_spectrum(energies, roots)

            if self.opts.verbose:
                sys.stderr.write("{0}: Initialized spectrum [{1: .3f} : {2: .3f}] ({3} points)\n".format(pname, energies[0], energies[-1], len(energies)))

            if self.opts.header:
                self.write("{0} Roots were Lorentzian broadened with width {1} eV ".format(self.opts.cchar, self.opts.width))
                self.write("{0} Spectrum generated with {1} points.".format(self.opts.cchar, self.opts.npoints))

        else:
            spectrum = roots
            if self.opts.header:
                self.write("{0} Roots not broadened".format(self.opts.cchar))
                self.write("{0} No spectrum generated: output is list of raw excitations".format(self.opts.cchar))


        ## finally, dump data to stdout
        self.dump_data(spectrum)

        if self.opts.verbose:
            sys.stderr.write ("{0}: Done\n".format(pname))
        


def main():
    ##
    ## Parse command line options.  Note we later make all lowercase,
    ## so from the user's point of view they are case-insensitive.
    ##
    usage = "%prog [options]\n\n"
    desc = "Reads NWChem output from stdin or specified infile, parses for the linear response TDDFT or DFT vspec excitations, and prints the absorption spectrum to stdout or writes it to outfile.  It will optionally broaden peaks using a Lorentzian with FWHM of at least two energy/wavelength spacings.  By default, it will automatically determine data format (tddft or vspec) and generate a broadened spectrum in eV."
    desc_wrap = textwrap.wrap(desc, 80)
    for s in desc_wrap:
        usage += s + "\n"

    example = 'Create absorption spectrum in nm named "spectrum.dat" from the NWChem output file "water.nwo" named spectrum.dat with peaks broadened by 0.3 eV and 5000 points in the spectrum.'
    example_wrap = textwrap.wrap(example, 80)
    usage += "\nExample:\n\n\t"+"nw_spectrum -b0.3 -p5000 -wnm < water.nwo > spectrum.dat\n\n\tor:\n\n\tnw_spectrum -b0.3 -p5000 -wnm -i water.nwo -o spectrum.dat\n\n"
    for s in example_wrap:
        usage += s + "\n"

    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--infile", type="string", dest="infile",
                      help="input file (reads from stdin if unspecified)", metavar="FILE")
    parser.add_option("-o", "--outfile", type="string", dest="outfile",
                      help="output file (writes to stdout if unspecified)", metavar="FILE")
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


    parser.set_defaults(width=0.1, npoints=2000, units="eV", datafmt="auto",
                        delim="    ", makespec=True, header=True, cchar="#",
                        nbin=20, verbose=False)
    (opts, args) = parser.parse_args()

    if opts.infile:
        instream = open(opts.infile, "r")
    else:
        instream = sys.stdin

    if opts.outfile:
        outstream = open(opts.outfile, "w")
    else:
        outstream = sys.stdout
        
    S = SpectrumProcessor(opts, instream, outstream)
    S.process_all()
    



if __name__ == "__main__":
    main()


