#!/usr/bin/python

# Amila Weerasinghe and Matthew Wyczalkowski
# amila@wustl.edu

import argparse
import pysam
from pybedtools import BedTool

def main():
    # from optparse import OptionParser
    usage_text = """usage: %(prog)s [options] BAM 
    takes a bam file and an adequately coverage cutoff to generate 
    a bed file containing the adequately covered regions"""

    parser = argparse.ArgumentParser( description=usage_text )
    parser.add_argument( '--version', action='version', version='%(prog)s 0.1' )
    parser.add_argument( "-d", dest="delta", default="15", help="minimum adequate read depth" )
    parser.add_argument( "-o", dest="outFile", default="./sample_bed.bed", help="output bed file name" )
    parser.add_argument( "bam", help="bam filename" )
    # parser.add_argument( "--bam", dest="bamFile", default='sample.bam', help="bam filename" )

    args = parser.parse_args()
    # print( args, args.bam, len(args.bam) )
    # if ( len(args.bam) != 1 ):
    #     parser.error( "Exactly one BAM file should be provided." )
    
    generateBedfileFromBam( args.outFile, args.bam, int(args.delta) )

def generateBedfileFromBam( outFile, bamFile, delta ):
    # for the provided BAM file, generate bed files indicating adequately covered regions
    # print( "bam files given=", args.bamFiles )
    ##TODO check whether bams are indexed and sorted, then perform indexing or sorting if needed
    print( "generating adequately covered bed file for ", bamFile )
    bedFile = BedTool( bamFile ).genome_coverage( bg=True ) # bg=True gives read depth in bed graph format
        
    filteredBedFile = bedFile.filter( lambda x: int(x.name) >= delta ) # bedfile object carries the fourthcolumn as attribute "name"  # keep only the regions with read depth >= adequate sample coverage "delta" 
    # filteredBedFile.saveas( args.outFile )
        
    mergedBedFile = filteredBedFile.merge() # merge adequately covered regions together
    mergedBedFile.saveas( outFile )

    print( "bed file generation is done, written to ", outFile )
    

if __name__ == '__main__':
    main()