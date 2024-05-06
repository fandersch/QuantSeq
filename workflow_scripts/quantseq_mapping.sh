#!/bin/bash
IN="$1"
REF="$2"
OUT="$3"
ThreeUTR_BED="$4"
LOGS="logs/quantseq_mapping"

rm -r $OUT
rm -r ${LOGS}

mkdir -p $OUT
mkdir -p $OUT/map
mkdir -p $OUT/filter
mkdir -p $OUT/summary
mkdir -p ${LOGS}

jobName='quantseq_mapping'
for file in ${IN}/*.fq.gz;
do
    echo $file
    fname=`basename ${file%.fq.gz}`
    fname_output=${fname}.fq_slamdunk_mapped.bam
    echo $fname_output

    # create the script for each sample
    script=${LOGS}/${fname}_${jobName}.sh
    cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --time=120
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o ${LOGS}/$fname.out
#SBATCH -e ${LOGS}/$fname.err
#SBATCH --job-name $jobName

slamdunk map -r $REF -o $OUT/map -5 12 -n 100 -t 16 -q -ss $file
slamdunk filter -o $OUT/filter -b $ThreeUTR_BED -t 16 "$OUT/map/${fname_output}"
alleyoop summary -o $OUT/summary/${fname}.txt $OUT/filter/${fname}*.bam
EOF

    cat $script;
    sbatch $script
done
