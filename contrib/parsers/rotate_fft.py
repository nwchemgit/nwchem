#!/usr/bin/python

import sys
import math


def rotate_spectrum (data):
    for v in data:
        w = v[0]
        re = v[1]
        im = v[2]
        ab = v[3]
        
        r = math.sqrt (re**2 + im**2)
        if abs (r - ab) > 1e-6:
            raise Exception ("abs not equal to sqrt(re^2 + im^2)")

        angle = abs (math.atan2 (im, re))

        if angle > math.pi:
            raise Exception ("atan2 out of range")

        re_out = ab * math.cos (angle)
        im_out = ab * math.sin (angle)

        yield [w, re_out, im_out, ab]

        
def parse_stdin ():
    lines = sys.stdin.readlines ()
 
    # strip comments: ignore any line with "#" anywhere in it
    rawinp = [ l.strip().split() for l in lines if not "#" in l ]

    # convert data to floats
    data = [ [ float(v) for v in l ] for l in rawinp]
    
    return data


def main ():
    data = parse_stdin ()
    data_rot = rotate_spectrum (data)

    for d in data_rot:
        print ("%20.10e\t%20.10e\t%20.10e\t%20.10e" %(d[0], d[1], d[2], d[3]))


if __name__ == "__main__":
    main()
