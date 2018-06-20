#!/usr/bin/python

# Amila Weerasinghe and Matthew Wyczalkowski
# amila@wustl.edu

import argparse
import pysam
from pybedtools import BedTool
import numpy as np

def main():
    # from optparse import OptionParser
    usage_text = """usage: %(prog)s [options] bedFile1 bedFile2 ... bedFileN ...
    takes a set of bed files and a population adequately coverage cutoff to generate 
    a bed file containing the adequately covered regions in the population"""

    parser = argparse.ArgumentParser( description=usage_text )
    parser.add_argument( '--version', action='version', version='%(prog)s 0.1' )
    parser.add_argument( "-m", dest="m_ratio", default="0.75", help="cutoff for the adequately covered sample ratio" )
    parser.add_argument( "-o", dest="outFile", default="stdout", help="output bed file name" )
    parser.add_argument( "bed", nargs='+', help="bed filenames" )
    # parser.add_argument( "--bed", dest="bedFiles", metavar='sample1.bed', type=str, nargs='+', help="bed filenames" )

    args = parser.parse_args()

    if ( len(args.bed) == 0 ):
        parser.error( "At least one bed file should be provided." )
    bedFiles = args.bed
    # print( args )
    print( "bed files given=", bedFiles )
    # numSamples = len( bedFiles )
    getAdequatelyCoveredRegionsInPopulation( bedFiles, args.m_ratio, args.outFile )
    

def getAdequatelyCoveredRegionsInPopulation( bedFiles, m_ratio, outFile ):
    intervalInfo = generateArraysOfStartEndPositions( bedFiles ) # a dictionary to store interval-array for chromosome number keys
    segmentsBed = createBedFileFromSegments( intervalInfo ) # a BedTool object containing segments in bed format
    # print( "segmentsBed= ", segmentsBed )
    popCovBed = getPopulationCoverage( segmentsBed, bedFiles ) # a BedTool object containing segments in bed format annotated with the number of samples covering those segments
    # print( popCovBed )
    cutoffNumberOfSamples = float(m_ratio)*len(bedFiles)
    # print( cutoffNumberOfSamples )
    filteredPopCovBed = popCovBed.filter( lambda x: int(x.name) >= cutoffNumberOfSamples ) # bedfile object carries the fourthcolumn as attribute "name"  # keep only the regions with num_samples/total num_samples >= ratio 
    # print( filteredPopCovBed )
    filteredPopCovBed = filteredPopCovBed.merge() # merge adjacent regions to be more cleaner

    if outFile == "stdout":
        print( filteredPopCovBed )
    else:
        filteredPopCovBed.saveas( outFile )

def generateArraysOfStartEndPositions( bedFiles ):
    intervalInfo = {}
    for bedFile in bedFiles:
        bed = BedTool(bedFile)
        for feature in bed:
            #print(feature)
            # print(feature.chrom, feature.start, feature.stop)
            # add the start and end values to the array
            if feature.chrom not in intervalInfo.keys():
                intervalInfo[feature.chrom] = list()
            intervalInfo[feature.chrom].append( feature.start )
            intervalInfo[feature.chrom].append( feature.stop-1 )
    # now take each array for each chrom, take unique and sort the positions            
    for chrom in intervalInfo.keys():
        intervalInfo[chrom] = np.unique( np.sort( np.array( intervalInfo[chrom] ) ) )
    return intervalInfo

def createBedFileFromSegments( intervalInfo ):
    segments = [] # an empty string to store info in bed-format
    for chrom in intervalInfo.keys():
        starts = intervalInfo[chrom][:-1]
        stops = intervalInfo[chrom][1:]
        for start, stop in zip( starts, stops ): # zip makes an iterator object. Faster than making an intermediate array
            segments.append( '\t'.join( (str(chrom), str(start), str(stop+1), "0" ) ) )# stop+1 b/c end is 1-base; added a zero to the end to store number of samples downstream
    # create the BedTool object using the segments string
    segmentsBed = BedTool( '\n'.join( segments ), from_string=True )
    return segmentsBed

def checkIfFeatureIsCoveredInSample( feature, sampleBed ):
    # print( "Testing feature= ", feature, "feature name= ", feature.name )
    if sampleBed.any_hits( feature ):
        # if feature.name=='0': # this is the first time this segment is considered
        #     feature.name = '1'
        # else:
        feature.name = str( int(feature.name) + 1 ) # name attribute is the fourth column
        # print( "\t There's a hit. Name increased by one " )
    # else:
    #     if feature.name=='.': # this is the first time this segment is considered (need to set to zero bc there can be regions not covered by any sample)
    #         feature.name = '0'
    #     print( "\t No hit. Name kept at ", feature.name )
    # print( "\t New feature= ", feature, "feature name= ", feature.name )
    return feature  

def getPopulationCoverage( segmentsBed, bedFiles ): 
    for sampleBedFile in bedFiles:
        sampleBed = BedTool( sampleBedFile )
        # print( sampleBed )
        segmentsBed = segmentsBed.each( checkIfFeatureIsCoveredInSample, sampleBed ) 
        # for feature in segmentsBed:
        #     newFeature = checkIfFeatureIsCoveredInSample( feature, bed )
        #     print( newFeature )
            # newFeatureBed = BedTool( newFeature )
            # segmentsBedNew.cat( newFeatureBed )
        
        # segmentsBed = segmentsBedNew
    return segmentsBed

            
if __name__ == '__main__':
    main()