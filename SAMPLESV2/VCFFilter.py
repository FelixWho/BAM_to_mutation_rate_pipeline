import argparse
import vcf

def main():
    usage_text = """usage: %(prog)s [options] VCF Filter ...
    Filters a VCF file by mutation type, returning a new VCF File"""

    parser = argparse.ArgumentParser( description=usage_text )
    parser.add_argument( '--version', action='version', version='%(prog)s 0.1' )
    parser.add_argument( dest='vcf', metavar='VCF_file', type=str, help='VCF file to filter' )
    parser.add_argument( "-o", default = "filtered.vcf", dest="outFile", help="output VCF file name" )
    parser.add_argument( "-f", dest = 'filter', default = "ALL", help = "filter argument e.g. \'SNP\'" )

    args = parser.parse_args()
    filterVCF(args.filter.upper(), args.vcf, args.outFile)

def filterVCF(filterBy, VCFFile, outFile):
    vcf_reader = vcf.Reader(open(VCFFile))
    vcf_writer = vcf.Writer(open(outFile, 'w'), vcf_reader)

    for record in vcf_reader:
        if not hasattr(record, "INFO"):
            print("INFO not found for record at position: %d" % record.POS)
            continue
        if filterBy == "ALL":
            vcf_writer.write_record(record)
            continue
        if record.INFO["VT"][0] == filterBy:
                vcf_writer.write_record(record)

if __name__ == '__main__':
    main()