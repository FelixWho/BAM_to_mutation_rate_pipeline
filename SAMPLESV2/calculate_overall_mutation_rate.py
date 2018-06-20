#!/usr/bin/python

# Felix Hu Attempt 1

import argparse
import pysam
from pybedtools import BedTool
import numpy as np
import vcf

def main():
    usage_text = """usage: %(prog)s [options] PopulationCoveredROI VCF ...
    takes a population covered ROI file and a vcf file containing the mutation calls 
    from all the samples in the population to calculate the mutation rate for all ROIs"""

    parser = argparse.ArgumentParser( description=usage_text )
    parser.add_argument( '--version', action='version', version='%(prog)s 0.1' )
    parser.add_argument( "-o", dest="outFile", default="stdout", help="output file name" )
    parser.add_argument( dest='ROI', metavar='ROI_file', type=str, help='population covered ROI file name' )
    parser.add_argument( dest='vcf', metavar='VCF_file', type=str, help='VCF file containing mutation calls for all the samples in the population' )
    # parser.add_argument( "--bed", dest="bedFiles", metavar='sample1.bed', type=str, nargs='+', help="bed filenames" )

    args = parser.parse_args()

    ROIBed = BedTool( args.ROI ) 
    
    result = intersectMutationRegardlessOfMutationType( args.vcf )

    mutation_count = intersectBedsAndRetNumMutations(result[0], ROIBed, result[1])
    num_samples = result[2]
    track_length = ROIBed.sort().total_coverage()
    print( "\n\n"+str(float(mutation_count / (num_samples * track_length))) + "\n\n")

def intersectMutationRegardlessOfMutationType(vcfFile):
    vcfData = vcf.Reader( open(vcfFile, 'r') )

    mut_data = dict()
    bed_data = list()
    num_samples_found = False
    num_samples = -1
    for record in vcfData:
        #print(str(record.num_called))
        if not num_samples_found:
            num_samples = len(record.samples)
            num_samples_found = True
        featureString = '\t'.join( [str(record.CHROM), str(record.start), str(record.end)] )
        mut_data[(record.start, record.end)] = num_samples - len(record.get_hom_refs())
        bed_data.append(featureString)

    return ( BedTool('\n'.join(bed_data), from_string=True), mut_data, num_samples )

def intersectBedsAndRetNumMutations(vcfBed, ROIBed, mut_data):
    intersectedBed = ROIBed.intersect(vcfBed)
    total_mutations = 0
    for interval in intersectedBed:
        total_mutations = total_mutations + mut_data[(interval.start, interval.stop)]

    return total_mutations


if __name__ == '__main__':
    main()