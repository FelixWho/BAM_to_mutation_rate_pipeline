import vcf
import argparse
import copy
import collections
from random import uniform
from random import randint

def main():
    # from optparse import OptionParser
    usage_text = """usage: %(prog)s [options] MUT
    takes a vcf file to generate a vcf file 
    containing randomly induced mutations"""

    parser = argparse.ArgumentParser( description=usage_text )
    parser.add_argument( '--version', action='version', version='%(prog)s 0.1' )
    parser.add_argument( "-d", dest="delta", default="20.0", help="probability at any single non-mutated position of inducing a mutation as a percentage" )
    parser.add_argument( "-o", dest="outFile", default="./sample_vcf.vcf", help="output vcf file name" )
    parser.add_argument( "v", help="vcf filename" )
    
    args = parser.parse_args()

    induce_mutations(args.v, args.outFile, float(args.delta)/10.0)
    # print( args, args.bam, len(args.bam) )
    # if ( len(args.bam) != 1 ):
    #     parser.error( "Exactly one vcf file should be provided." )

def induce_mutations (inFile, outFile, delta):
    
    vcf_reader = vcf.Reader(open(inFile))
    vcf_writer = vcf.Writer(open(outFile, 'w'), vcf_reader)

    for record in vcf_reader:
        rec_toWrite = copy.deepcopy(record)
        print(record.num_called)
        if record.FORMAT.split(":")[0] != "GT":
            print("Error with FORMAT column at POSITION=" + str(record.POS) + ": \'GT\' is missing.\n")
            continue

        for sm in range(len(record.samples)):
            gt_read = str(record.genotype(record.samples[sm].sample)["GT"])
            f_vals = [record.samples[sm].data[vx] for vx in range(len(record.FORMAT.split(":")))]
            f_keys = record.FORMAT.split(":")
            rec_toWrite.samples[sm].data = collections.namedtuple('CallData', f_keys)
            if not already_mutated(gt_read):
                if uniform(0,9) < delta:
                    mutation_type = randint(0, 2)
                    mut_type_str = ""
                    
                    if mutation_type == 0:
                        mut_type_str = str(randint(1,len(record.ALT))) + ("|" if gt_read[1] == "|" else "/") + "0"
                    elif mutation_type == 1:
                        mut_type_str = "0" + ("|" if gt_read[1] == "|" else "/") + str(randint(1,len(record.ALT)))
                    else:
                        mut_type_str = str(randint(1,len(record.ALT))) + ("|" if gt_read[1] == "|" else "/") + str(randint(1,len(record.ALT)))

                    f_vals[0] = mut_type_str
            rec_toWrite.samples[sm].data = rec_toWrite.samples[sm].data._make(f_vals)
        
        vcf_writer.write_record(rec_toWrite)

def already_mutated(gt):
    if len(gt) != 3:
        return False

    if gt[1] == "/":
        if gt == "0/0":
            return False
    else:
        if gt == "0|0":
            return False
    return True

if __name__ == '__main__':
    main()