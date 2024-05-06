#!/bin/bash
IN="$1"
OUT="$2"
ThreeUTR_GTF="$3"
LOGS="logs/featureCounts"

rm -r $OUT
rm -r ${LOGS}

mkdir -p $OUT
mkdir -p ${LOGS}

jobName='featureCounts'
outdir=`basename ${OUT}`

# create the script for each sample
script=${LOGS}/${jobName}_${outdir}.sh
cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=60
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o ${LOGS}/${jobName}_${outdir}.out
#SBATCH -e ${LOGS}/${jobName}_${outdir}.err
#SBATCH --job-name $jobName

featureCounts -s 1 -M -T 1 -a $ThreeUTR_GTF -o ${OUT}/slamdunk_3UTRs.count $IN/*.bam
EOF

cat $script;
sbatch $script
