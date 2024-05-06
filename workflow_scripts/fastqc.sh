#!/bin/bash
DIR=`pwd`
INPUT_PATH="${DIR}/$1"
OUTPUT_DIR="${DIR}/$2"
LOGS="${DIR}/logs/fastqc"

rm -r "$OUTPUT_DIR/fastqc"
rm -r ${LOGS}
mkdir -p "$OUTPUT_DIR/fastqc"
mkdir -p ${LOGS}

jobName='fastqc'
jobids=()

for file in ${INPUT_PATH}/*.*;
do
    echo $file
    fname=`basename ${file%.bam}`
    fname=${fname%.fastq.gz}
    fname=${fname%.fq.gz}
    fname=${fname%.fastq}
    fname=${fname%.fq}
    
    # create the script for each sample
    script_fastqc=${LOGS}/${fname}_${jobName}.sh
    cat <<EOF > $script_fastqc
#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1     # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --mem=30G  # memory in Mb
#SBATCH --time=30  # time
#SBATCH -o ${LOGS}/${fname}.out # STDOUT
#SBATCH -e ${LOGS}/${fname}.err # STDERR
#SBATCH --job-name $jobName

fastqc ${file} -o "${OUTPUT_DIR}/fastqc"

EOF
    cat $script_fastqc
    jobid_buff=$(sbatch $script_fastqc)
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
script_multiqc=${LOGS}/${jobName}.sh
cat <<EOF > $script_multiqc
#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1     # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --mem=10G  # memory in Mb
#SBATCH --time=30  # time
#SBATCH -o ${LOGS}/$jobName.out # STDOUT
#SBATCH -e ${LOGS}/$jobName.err # STDERR
#SBATCH --job-name $jobName

cd ${OUTPUT_DIR}
multiqc ${OUTPUT_DIR}/fastqc

EOF

cat $script_multiqc;
sbatch --dependency=afterok:${jobids_str} $script_multiqc
