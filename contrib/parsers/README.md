### Script for extracting data from RT-TDDFT simulations

The following options can be used to extract data from RT-TDDFT simulations using the `nw_rtparse.py` script
```
$ python nw_rtparse.py  -help
Usage: nw_rtparse [options] output.nwo [output2.nwo] 

Parse NWChem real-time TDDFT output for time-dependent quantities.


Options:
  -h, --help            show this help message and exit
  -t STR, --tag=STR     parse for simulation tag STR, defaults to '<rt_tddft>'
  -g STR, --geometry=STR
                        extract data for geometry STR, defaults to 'system'
  -x STR, --extract=STR
                        target data to extract: dipole (default), efield,
                        energy, S2, charge
  -p STR, --polarization=STR
                        target polarization: x (default), y, z
  -s STR, --spin=STR    target spin: closedshell (default), alpha, beta, total
  -C, --clean           clean output; data only, no header or comments
  -d STR, --delim=STR   use STR as output separator (four spaces default)
  -z, --zero            zero data at t=0
  -R VAL, --tolerance=VAL
                        tolerance for checks (default 1e-5)
  -c, --compare         read in two files and compare
  ```
