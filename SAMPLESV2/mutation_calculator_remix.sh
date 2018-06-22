#!/bin/bash

# USAGE

# $ bash mutation_calculator_remix.sh [-d delta_bam_to_bed] [-D delta_bed_to_adequate_covered_regions] [-f vcf_filter_parameter] -v vcf_file -r roi_file [-R another_roi_file] [-o output_file] BAM_1 [BAM_2 ... ]
  # -d specifies read depth of BAM file required for DNA sequence to be included in BED file                            (optional)
  # -D specifies depth of coverage required for BED interval to be included in population-adequate-coverage BED file    (optional)
  # -f specifies what mutation type to filter for in VCF file i.e. "SNP" "INDEL" "ALL"                                  (optional)
  # -v specifies path of VCF file                                                                                       (required)
  # -r specifies path of ROI (region-of-interest) file                                                                  (required)
  # -R specifies an optional ROI file for comparison with first ROI file                                                (optional)
  # -o specifies path of output file; default is standard output                                                        (optional)
  # BAM1 specifies path of a required BAM file                                                                          (required)
  # BAM2 ... specifies other BED files                                                                                  (optional)

num_beds=0
delta_bam_to_bed=15
delta_bed_to_adequate_covered_regions=0.75
bam_inputted=false
vcf_inputted=false
roi_inputted=false
optional_roi_inputted=false
result_to_stdout=true
filter=all

while getopts ":d:D:f:v:r:R:o:" opt; do
  case $opt in
    d)
      echo "delta_bam_to_bed set to: $OPTARG" >&2
      delta_bam_to_bed=$OPTARG
      ;;
    D)
      echo "delta_bed_to_adequate_covered_regions set to: $OPTARG" >&2
      delta_bed_to_adequate_covered_regions=$OPTARG
      ;;
    f)
      echo "filtering VCF by: $OPTARG" >&2
      vcf_filter_parameter=$OPTARG
      ;;
    v)
      echo "taking in VCF file: $OPTARG" >&2
      vcf_file=$OPTARG
      if $vcf_inputted
      then
        vcf_inputted=false
      else
        vcf_inputted=true
      fi
      ;;
    r)
      echo "taking in ROI file: $OPTARG" >&2
      roi_file=$OPTARG
      roi_inputted=true
      ;;
    R)
      echo "taking in optional ROI file: $OPTARG" >&2
      optional_roi_file=$OPTARG
      optional_roi_inputted=true
      ;;
    o)
      echo "result will go to file: $OPTARG"
      output_file=$OPTARG
      result_to_stdout=false
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))

if ! $vcf_inputted
then
    echo "Must input a single VCF file"
    exit 1
fi
if ! $roi_inputted
then
    echo "Must input a ROI file"
    exit 1
fi
if ! (($#>0));
then
    echo "Must input at least one BAM file"
    exit 1
fi

in_files=""
for ((i=1;i<=$#;i++)); do
  in_files=$in_files" bed_"$i".bed"
  out_file="bed_"$i".bed"
  python ./bam_to_bed.py "${!i}" -o $out_file -d $delta_bam_to_bed
done

out_pac_file="adequate_covered_regions.bed"
intersected_file_1="intersect_roi_and_pac.bed"
intersected_file_2="intersect_roi_and_pac_2.bed"
filtered_vcf="filtered.vcf"

#filter VCF file
python ./VCFFilter.py $vcf_file  -o $filtered_vcf -f $vcf_filter_parameter

#create adequately covered regions bed
python ./beds_to_covered_regions.py $in_files -o $out_pac_file -m $delta_bed_to_adequate_covered_regions

#intersect adequately covered regions bed and Regions Of Interest bed(s)
python ./intersect_roi_and_popCovBed.py -o $intersected_file_1 $roi_file $out_pac_file

if $optional_roi_inputted
then
  python ./intersect_roi_and_popCovBed.py -o $intersected_file_2 $optional_roi_file $out_pac_file
fi

#calculate mutation rate(s)
RESULT_1=$(python ./calculate_overall_mutation_rate.py $intersected_file_1 $filtered_vcf 2>&1)
if $result_to_stdout
then
  echo -n "$roi_file has mutation rate per base pair: "
  echo "scale=10 ; $RESULT_1 / 1" | bc
else
  echo -n "$roi_file has mutation rate per base pair: " >> $output_file
  echo "scale=10 ; $RESULT_1 / 1" | bc >> $output_file
fi

if $optional_roi_inputted
then
  RESULT_2=$(python ./calculate_overall_mutation_rate.py $intersected_file_2 $filtered_vcf 2>&1)
  if $result_to_stdout
  then
    echo -n "$optional_roi_file has mutation rate per base pair: "
    echo "scale=10 ; $RESULT_2 / 1" | bc
    echo -n "mutation rate ratio between $roi_file and $optional_roi_file: "
    echo "scale=10 ; $RESULT_1 / $RESULT_2" | bc
  else
    echo -n "$optional_roi_file has mutation rate per base pair: " >> $output_file
    echo "scale=10 ; $RESULT_2 / 1" | bc >> $output_file
    echo -n "mutation rate ratio between $roi_file and $optional_roi_file: " >> $output_file
    echo "scale=10 ; $RESULT_1 / $RESULT_2" | bc >> $output_file
  fi
fi
