#!/usr/bin/python

# Amila Weerasinghe and Matthew Wyczalkowski
# amila@wustl.edu

import argparse
import pysam
from pybedtools import BedTool
import numpy as np

def main():
    # from optparse import OptionParser
    usage_text = """usage: %(prog)s [options] ROI_file PopulationCoverage_file ...
    takes an ROI file and a bed file containing the adequately covered regions in the population
    and return the intersect of those two in bedfile format"""

    parser = argparse.ArgumentParser( description=usage_text )
    parser.add_argument( '--version', action='version', version='%(prog)s 0.1' )
    parser.add_argument( "-o", dest="outFile", default="stdout", help="output bed file name" )
    parser.add_argument( dest='ROI', metavar='ROI_file', type=str, help='ROI file name' )
    parser.add_argument( dest='popCovBed', metavar='Bed_file', type=str, help='Bed file containing adequately covered regions in the population' )
    # parser.add_argument( "--bed", dest="bedFiles", metavar='sample1.bed', type=str, nargs='+', help="bed filenames" )

    args = parser.parse_args()
    print( args )

    ROIBed = BedTool( args.ROI ) 
    popCovBed = BedTool( args.popCovBed )

    intersectedBed = ROIBed.intersect( popCovBed )

    if args.outFile == "stdout":
        print( intersectedBed )
    else:
        intersectedBed.saveas( args.outFile )
    

            
if __name__ == '__main__':
    main()