#!/bin/bash
#######################################################################################################################################################################################################################
# script to trim adaptors for fastq files and run fastqc after trimming
# two modules possible: cutadapt which need the adaptor sequences to be speficied
# but the adavantage of using cutadapt isthat it can trim polyA at the same time (however this polyA trimming is not
# implemented for some reasons)
# trimglore which does not need the adaptor sequences and this option is removed now
#
#######################################################################################################################################################################################################################

#adaptor_seq="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG" # quant-seq adaptor
#adaptor_seq="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"; # Nextera adatptor
#trim_polyA="TRUE";

DIR=`pwd`
DIR_input="${DIR}/$1"
DIR_trimmed="${DIR}/$2"
DIR_fastqc="${DIR}/$2/FASTQCs_trimmed/"
base_quality_threshold=$3
firstbpToClip=$4
adapter_stringency=$5
LOGS="$DIR/logs/trimming_fastq/"

rm -r $DIR_trimmed
rm -r $DIR_fastqc
rm -r "${DIR}/logs/trimming_fastq"

mkdir -p $DIR_trimmed
mkdir -p $DIR_fastqc
mkdir -p "${DIR}/logs/trimming_fastq"

jobName='trimming'
for file in $DIR_input/*;
do
    echo $file
    fname=`basename ${file%.fq.gz}`
    fname=${fname%.fastq.gz}
    fname=`basename ${fname%.fq}`
    fname=${fname%.fastq}
    fname_input=${fname}

    IFS='#' read -ra SAMPLE_ID <<< "${fname}"
    if [ ${#SAMPLE_ID[@]} -eq 2 ]; then
        echo "${SAMPLE_ID[1]}"
        fname=${SAMPLE_ID[1]}
    else
        echo "${SAMPLE_ID[0]}"
        fname=${SAMPLE_ID[0]}
    fi


    cd $DIR/logs/trimming_fastq || exit
    out="${fname%.R1}_trimmed.fq.gz";
    echo ${out}

    # creat the script for each sample
    script=$PWD/${fname}_${jobName}.sh
    cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --time=120
#SBATCH --mem=15G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $LOGS/$fname.out
#SBATCH -e $LOGS/$fname.err
#SBATCH --job-name $jobName

if [ "${firstbpToClip}" -eq "0" ]; then
  trim_galore $file -j 8 --quality ${base_quality_threshold} --fastqc_args "--outdir ${DIR_fastqc}" --output_dir ${DIR_trimmed} --stringency $adapter_stringency
else
  trim_galore $file -j 8 --quality ${base_quality_threshold} --clip_R1 ${firstbpToClip} --fastqc_args "--outdir ${DIR_fastqc}" --output_dir ${DIR_trimmed} --stringency $adapter_stringency
fi

if [ "${DIR_trimmed}/${fname_input}_trimmed.fq.gz" != "${DIR_trimmed}/${out}" ]; then
  mv ${DIR_trimmed}/${fname_input}_trimmed.fq.gz ${DIR_trimmed}/${out}
fi

EOF
    cat $script
    jobid_buff=$(sbatch $script)
    if [[ "$jobid_buff" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
      jobid_buff="${BASH_REMATCH[1]}"
      jobids+=("$jobid_buff")
    else
      echo "sbatch failed"
    fi
done

jobids_str=$( IFS=$':'; echo "${jobids[*]}" )
echo $jobids_str

jobName='multiqc'
# create the script for each sample
script_multiqc=$LOGS/${jobName}.sh
cat <<EOF > $script_multiqc
#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1     # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --mem=10G  # memory in Mb
#SBATCH --time=60
#SBATCH -o $LOGS/$jobName.out # STDOUT
#SBATCH -e $LOGS/$jobName.err # STDERR
#SBATCH --job-name $jobName

cd ${DIR_trimmed}
multiqc ${DIR_fastqc}

EOF

cat $script_multiqc;
sbatch --dependency=afterok:${jobids_str} $script_multiqc
