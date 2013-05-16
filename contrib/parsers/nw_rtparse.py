##
## nw_rtparse
##
## Kenneth Lopata
## Last modified: 2013-05-01
##
## Python script for parsing NWChem real-time TDDFT output for various
## time-dependent quantities (e.g., dipole moment, electric field, etc)
##
## nw_rtparse --help
##
##

import sys
from optparse import OptionParser
from signal import signal, SIGPIPE, SIG_DFL


def check_version ():
    """Check version and ensure new print and string formatting is
    allowed.  Raises and exception if not satisfied."""

    if sys.version_info < (2, 6):
        raise Exception("This script requires python >= 2.6")

    try:
        newstring = "oldstring {v}".format(v=3.14)
    except:
        raise Exception("This script requires string.format()")


def parse_input (opts, labels, dataindex, lines):
    
    data = []
    for l in lines:

        if all ([tag in l for tag in labels]):

            try:
                vals = l.strip().split()
            except:
                raise Exception ("Failed to parse line: {0}".format(l))

            # tag t v1 v2 v3 ...
            try:
                tagin = str(vals[0])  
                timein = float(vals[1])
                valin = float(vals[dataindex])
            except:
                raise Exception ("Failed to convert data (data index {0}) from line: {1}".format(dataindex, l))
            
            data.append ([timein, valin])


    ndata = len(data)

    if ndata <= 0:
        raise Exception ("Failed to find any data, check target and options.")
    
    if opts.header:
        print ("# Number of data points: {0}".format(ndata))
        print ("# Time range: [{0}: {1}]".format(data[0][0], data[-1][0]))
            
    return data



    

def postprocess_check (opts, data):
    
    tprev = -9999.0
    d0 = data[0][1]

    for d in data:
        t = d[0]
        
        if abs(t - tprev) < 1e-6:
            raise Exception ("Error: Repeating times found--ambiguous target.  Multiple geometries?  Using closedshell target for openshell or spin-orbit?")

        tprev = t

        if opts.zero:
            d[1] = d[1] - d0

    if opts.header and opts.zero:
#        print ("# [postprocess] zeroed at first time point")
        print ("# Postprocessing: zeroed at first time point")

    return data


def check_args_determine_labels (opts):
    """Check arguments and generate tags for parsing the output"""
    labels = []
    
    ## do these outputs have x,y,z values, geom, spin tags etc?
    vector = False
    spin = False
    geom = False
    tag = True  #all TD properties are tagged with "tag"
    
    ## data type
    if opts.target == "dipole":
        labels.append ("Dipole moment")
        vector = True
        spin = True
        geom = True
        
    elif opts.target == "efield":
        labels.append ("Applied E-field")
        vector = True
        spin = True
        geom = True

    elif opts.target == "energy":
        labels.append ("Etot")
        vector = False
        spin = False
        geom = False

    elif opts.target == "S2":
        labels.append ("S^2")
        vector = False
        spin = False
        geom = False
        
    elif opts.target == "charge":
        labels.append ("Charge")
        vector = False
        spin = True
        geom = True

    else:
        raise Exception ("Invalid target: {0}".format(opts.target))

    
    ## tag
    if tag:
        labels.append (opts.tag)
    

    ## polarization
    if vector:
        if opts.polarization == "x":
            dataindex = 2
        elif opts.polarization == "y":
            dataindex = 3
        elif opts.polarization == "z":
            dataindex = 4
        else:
            raise Exception ("Invalid polarization: {0}".format(polarization))
    else:
        dataindex = 2


    ## spin type
    if spin:
        if opts.spin == "closedshell":
            pass
        elif opts.spin == "alpha":
            labels.append ("(alpha spin)")
        elif opts.spin == "beta":
            labels.append ("(beta spin)")
        elif opts.spin == "total":
            labels.append ("(total spin)")
        else:
            raise Exception ("Invalid spin type: {0}".format(opts.spin))
    else:
        pass

    
    ## geom
    if geom:
        labels.append (opts.geom)
    else:
        pass

        
    return labels, dataindex


def compare_data (opts, data1, data2):
    n = min (len(data1), len(data2))

    if opts.header:
        print ("")
        print ("Comparison")
        print ("----------")

    stat = True

    for i in range(n):
        t1 = data1[i][0]
        t2 = data2[i][0]

        v1 = data1[i][1]
        v2 = data2[i][1]

        tdiff = abs(t1-t2)
        vdiff = abs(v1-v2)

        if tdiff > opts.tolerance:
            print ("=> Problem: Inconsistent times: {0}, {1}".format(t1, t2))
            stat = False

        if vdiff > opts.tolerance:
            print ("=> Problem: Inconsistent values: {0}, {1}".format(v1, v2))
            stat = False


    if opts.header:
        print ("Compared {0} data points".format(n))
        print ("OK: {0}".format(stat))
        print ("")
            
    return stat


def print_output (opts, data):
    if opts.header:
#        outtag = "{0} {1} ({2})".format(opts.target, opts.polarization, opts.spin)
        outtag = "{0}".format(opts.target)
        
        print ("#")
        print ("#-----------------------------------------")
        print ("#   Time [au]        {0} [au]         ".format(outtag))
        print ("#-----------------------------------------")
        print ("#")
        
    for d in data:
        print ("{time: .10e}{delim}{val: .10e}".format(time=d[0], val=d[1], delim=opts.delim))



def main ():

    check_version ()


    # fixes broken pipe issue
    signal (SIGPIPE, SIG_DFL)
    
    ##
    ## parse command line args
    ##
    parser = OptionParser(usage="nw_rtparse [options] output.nwo [output2.nwo] \n\n"+
                          "Parse NWChem real-time TDDFT output for time-dependent quantities.\n")
    
    parser.add_option ("-t", "--tag", type="string", dest="tag",
                       help="parse for simulation tag STR, defaults to '<rt_tddft>'", metavar="STR")
    parser.add_option ("-g", "--geometry", type="string", dest="geom",
                       help="extract data for geometry STR, defaults to 'system'", metavar="STR")
    parser.add_option ("-x", "--extract", type="string", dest="target",
                       help="target data to extract: dipole (default), efield, energy, S2, charge", metavar="STR")
    parser.add_option ("-p", "--polarization", type="string", dest="polarization",
                       help="target polarization: x (default), y, z", metavar="STR")
    parser.add_option ("-s", "--spin", type="string", dest="spin",
                       help="target spin: closedshell (default), alpha, beta, total", metavar="STR")
    parser.add_option("-C", "--clean", action="store_false", dest="header",
                      help="clean output; data only, no header or comments")
    parser.add_option("-d", "--delim", type="string", dest="delim",
                      help="use STR as output separator (four spaces default)", metavar="STR")
    parser.add_option("-z", "--zero", action="store_true", dest="zero",
                      help="zero data at t=0")
    parser.add_option("-R", "--tolerance", type="float", dest="tolerance",
                      help="tolerance for checks (default 1e-5)", metavar="VAL")
    parser.add_option("-c", "--compare", action="store_true", dest="compare",
                      help="read in two files and compare")

    parser.set_defaults(tag = "<rt_tddft>", geom = "system", target = "dipole",
                        polarization="x", spin="closedshell", header = True,
                        delim = "    ", zero = False, tolerance=1e-5, compare = False)

    
    (opts, args) = parser.parse_args ()
    opts_dict = vars (opts)

    
     # check arguments and set up some parameters    
    labels, dataindex = check_args_determine_labels (opts)

    
    nargs = len(args)
    if opts.compare:
        if nargs != 2:
            raise Exception ("Exactly two arguments required if doing a comparison.")
    else:
        if nargs != 1:
            raise Exception ("Exactly one argument required for regular parse.")
        

    if opts.header:
        print ("# ======================================")
        print ("#  NWChem Real-time TDDFT output parser")
        print ("# ======================================")
        print ("# ")
        
        print ("# Runtime options")
        print ("# ---------------")
        for key, val in opts_dict.items ():
            print ("# {0} = {1}".format(key, val))

        

##
## read first data file
##
    fname1 = args[0]
    
    if opts.header:
        print ("#")
        print ("# File 1")
        print ("# ------")
        print ('# Filename: "{0}"'.format(fname1))

    try:
        file1 = open (fname1)
        lines1 = file1.readlines ()
        file1.close ()
    except:
        raise Exception ('Failed to read in data from file: '+'"'+fname1+'"')

    data1 = parse_input (opts, labels, dataindex, lines1)
    data1 = postprocess_check (opts, data1)

  
##
## read second data file (only for comparison) and carry out
## comparison
##  
    if opts.compare:
        fname2 = args[1]
        
        if opts.header:
            print ("#")
            print ("# File 2")
            print ("# ------")
            print ('# Filename: "{0}"'.format(fname2))

#        if fname1 == fname2:
#            raise Exception ("Invalid comparison between same files")

        try:
            file2 = open (fname2)
            lines2 = file2.readlines ()
            file2.close ()
        except:
            raise Exception ('Failed to read in data from file: '+'"'+fname2+'"')

        data2 = parse_input (opts, labels, dataindex, lines2)
        data2 = postprocess_check (opts, data2)
        stat = compare_data (opts, data1, data2)

        if stat:
            if opts.header:
                print ("Comparison passed ")
                print ("")
            else:
                print ('Comparison passed: "{0}" vs "{1}" with tolerance {3}'.
                       format(fname1, fname2, opts.target, opts.tolerance))
            sys.exit(0)
        else:
            if opts.header:
                print ("**********************")
                print ("* Comparison FAILED! *")
                print ("**********************")
                print ("")

            else:
                print ('Comparison FAILED: "{0}" vs "{1}" with tolerance {3}'.
                       format(fname1, fname2, opts.target, opts.tolerance))
            sys.exit(1)
            
        
    else:  #just a normal parse, print parsed data to stdout
        print_output (opts, data1)

        
if __name__ == "__main__":
    main()
