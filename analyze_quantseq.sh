#!/bin/bash
#SBATCH -N 1      # nodes requested
#SBATCH -n 1     # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --mem=20G  # memory in Mb
#SBATCH --job-name analyze_quantseq
mode=$1
umi=$2
gene_id=$3
rename=$4


if [[ ${mode} = "mm" ]]; then

  echo "mouse"
  UTRbed="/groups/zuber/USERS/florian.andersch/data/genomes/mouse/3UTR/mm10_refseq_ensembl_3UTR_fillup.bed"
  UTRgtf="/groups/zuber/USERS/florian.andersch/data/genomes/mouse/3UTR/mm10_refseq_ensembl_3UTR_fillup.gtf"
  genome="/groups/zuber/USERS/florian.andersch/data/genomes/mouse/fa/GRCm38.fa"

elif [[ ${mode} = "hs" ]]; then

  echo "human"
  UTRbed="/groups/zuber/USERS/florian.andersch/data/genomes/human/3UTR/hg38_refseq_ensembl_3UTR_fillup.bed"
  UTRgtf="/groups/zuber/USERS/florian.andersch/data/genomes/human/3UTR/hg38_refseq_ensembl_3UTR_fillup.gtf"
  genome="/groups/zuber/USERS/florian.andersch/data/genomes/human/fa/hg38_all.fa"

else
  echo "please specify mode: mm or hs"
  exit 1
fi

if [[ ! -z "$5" ]]; then
  UTRbed="$5"
fi

if [[ ! -z "$6" ]]; then
  UTRgtf="$6"
fi

if [[ ! -z "$7" ]]; then
  genome="$7"
fi

echo "UTRbed: ${UTRbed}"
echo "UTRgtf: ${UTRgtf}"
echo "Genome: ${UTRgtf}"


wait_for_jobs(){
  sleep 60  # seconds, give time to the schedular to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep florian. | grep -c $1) -eq 0 ]; then
      break
    else
      echo sleep another $((sleeptime / 60)) minutes...
      sleep $sleeptime
    fi
    n=$((n + 1))
    if [ $n -eq 5 ]; then
      sleeptime=300  # if still running after 10 mins, ask every 5 mins
    fi
    if [ $n -eq 10 ]; then
      sleeptime=600  # if still running after 30 mins, ask every 10 mins
    fi
  done
}

source ~/.bashrc
conda activate slamseq_quantseq

#######################
# rename fastq files  #
#######################
if [[ ${rename} = "rename" ]]; then
  srun --time=30 Rscript workflow_scripts/rename_filter_NGS_files.R fastq_${mode} input/input_${mode}.txt
fi

##############################################################################################################################################################################################################################################
# check read quality with fastq, @param1: directory of raw reads, @param2: input file wildcard to search for, @param3: output directory
##############################################################################################################################################################################################################################################

bash workflow_scripts/fastqc.sh fastq_${mode} QC_${mode}

output_dir=""
if [[ ${umi} = "umi" ]]; then

  ##############################################################################################################################################################################################################################################
  # extract umis from read and place into header, @param1: directory of merged reads, @param2: file wildcard (e.g. *.fq.gz), @param2: output directory
  ##############################################################################################################################################################################################################################################

  bash workflow_scripts/umitools_extract.sh fastq_${mode}/ fastq_umi_extracted_${mode}/
  wait_for_jobs umi_extr

  ##############################################################################################################################################################################################################################################
  # check read quality with fastq after umi extraction, @param1: directory of raw reads, @param2: input file wildcard to search for, @param3: output directory
  ##############################################################################################################################################################################################################################################

  bash workflow_scripts/fastqc.sh fastq_umi_extracted_${mode}/ QC_umi_${mode}

  output_dir="umi_extracted_"

fi

##############################################################################################################################################################################################################################################
# trimm adapter and first bases because of low quality,
# @param1: directory of reads (*.fq.gz),
# @param2: input file wildcard @param3: output directory,
# @param4: quality threshold,
# @param5: bases to trim at the beginning, @param6: adapter stringency
##############################################################################################################################################################################################################################################

bash workflow_scripts/trimming_fastq.sh fastq_${output_dir}${mode}/ trimmed_fastq_${output_dir}${mode}/ 25 4 3
wait_for_jobs trimming

##############################################################################################################################################################################################################################################
# map reads to genome with slam-dunk pipeline (NextGenMap), @param1: directory of trimmed umi reads (*.fq.gz), @param2: reference genome fats file, @param3: output directory, @param4: 3UTR bed-file
##############################################################################################################################################################################################################################################

#statements
bash workflow_scripts/quantseq_mapping.sh trimmed_fastq_${output_dir}${mode}/ ${genome} quantseq_mapping_${output_dir}${mode} ${UTRbed}

wait_for_jobs quantseq

if [[ ${umi} = "umi" ]]; then

  ##############################################################################################################################################################################################################################################
  # deduplicate reads that have the same sequence and the same UMI in the read name, @param1: directory of the mapped reads (*.bam), @param2: output directory
  ##############################################################################################################################################################################################################################################

  bash workflow_scripts/umitools_dedup.sh quantseq_mapping_umi_extracted_${mode}/filter/ quantseq_mapping_dedup_${mode}
  wait_for_jobs umi_dedu

  ##############################################################################################################################################################################################################################################
  # check read quality with fastq, @param1: directory of raw reads, @param2: input file type to search for, @param3: output directory
  ##############################################################################################################################################################################################################################################

  bash workflow_scripts/fastqc.sh quantseq_mapping_dedup_${mode}/ quantseq_mapping_dedup_${mode}/

  output_dir="dedup_"

fi

##############################################################################################################################################################################################################################################
#count reads per 3UTR gene region, @param1: directory of alignedreads (*.bam), @param2: output directory, param3: 3UTR gtf-file
##############################################################################################################################################################################################################################################
output_dir_feature_counts=""
subfolder="filter"
if [[ ${umi} = "umi" ]]; then
 output_dir_feature_counts="umi_"
 subfolder=""
fi

bash workflow_scripts/feature_counts_quantseq.sh quantseq_mapping_${output_dir}${mode}/${subfolder} featureCounts_${output_dir_feature_counts}${mode} ${UTRgtf}
wait_for_jobs featureC

##############################################################################################################################################################################################################################################
# multiqc
##############################################################################################################################################################################################################################################

bash workflow_scripts/final_multiQC_quantseq.sh . multiqc_final_${output_dir_feature_counts}${mode} ${mode}
wait_for_jobs final_mu

##############################################################################################################################################################################################################################################
# make tsv count file, @param1: input directory with slamdunk_3UTRs.count file in it
##############################################################################################################################################################################################################################################

bash workflow_scripts/gatherCounts_quantseq.sh featureCounts_${output_dir_feature_counts}${mode} ${gene_id} ${mode}
wait_for_jobs gatherCo

##############################################################################################################################################################################################################################################
# DESeq2 analysis
##############################################################################################################################################################################################################################################
srun --time=60 --mem=20G Rscript workflow_scripts/DESeq2_quantseq.R DEA_${mode} input/input_${mode}.txt featureCounts_${output_dir_feature_counts}${mode}/counts.tsv ${mode}
