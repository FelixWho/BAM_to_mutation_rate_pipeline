#!/bin/bash
echo "take bam, vcf, roi files . . . output mutation rate"

num_beds=0
delta_bam_to_bed=15
delta_bed_to_adequate_covered_regions=0.75
bam_inputted=false
vcf_inputted=false
roi_inputted=false
filter=all

while getopts ":d:D:f:b:v:r:" opt; do
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
      filter=$OPTARG
      ;;
    b)
      echo "taking in BAM file: $OPTARG" >&2
      bam_inputted=true
      ((num_beds++))
      out_file="bed"$num_beds".bed"
      python ./bam_to_bed.py $OPTARG -o $out_file -d $delta_bam_to_bed
      ;;
    v)
      echo "taking in VCF file: $OPTARG" >&2
      vcf_inputted=true
      vcf_file=$OPTARG
      ;;
    r)
      echo "taking in ROI file: $OPTARG" >&2
      roi_file=$OPTARG
      roi_inputted=true
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

if ! $bam_inputted
then
    echo "Must input at least 1 BAM file"
    exit 1
fi
if ! $vcf_inputted
then
    echo "Must input a VCF file"
    exit 1
fi
if ! $roi_inputted
then
    echo "Must input a ROI file"
    exit 1
fi

in_files=""
for ((i=1;i<=$num_beds;i++)); do
    in_files=$in_files" bed"$i".bed"
done

#create adequately covered regions bed
python ./beds_to_covered_regions.py $in_files -o "adequate_covered_regions.bed" -m $delta_bed_to_adequate_covered_regions

#intersect adequately covered regions bed and Regions Of Interest bed
python ./intersect_roi_and_popCovBed.py -o "intersect_roi_and_pac.bed" $roi_file "adequate_covered_regions.bed"

#filter VCF file
python ./VCFFilter.py $vcf_file  -o "filtered.vcf" -f $filter

#calculate mutation rate
python ./calculate_overall_mutation_rate.py "intersect_roi_and_pac.bed" "filtered.vcf"