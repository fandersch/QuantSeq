#!/bin/bash
DIR=`pwd`
IN="${DIR}/$1"
gene_id="$2"
mode="$3"
LOGS="${DIR}/logs/gatherCounts_quantseq"

rm -r ${LOGS}
mkdir -p ${LOGS}

jobName='gatherCounts_quantseq'

# create the script for each sample
script=${LOGS}/${jobName}_${1%/}.sh
cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=30
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o ${LOGS}/${jobName}_${1%/}.out
#SBATCH -e ${LOGS}/${jobName}_${1%/}.err
#SBATCH --job-name ${jobName}_${1%/}

Rscript /groups/zuber/USERS/florian.andersch/data/workspace/NGS_analysis_scripts/gatherCounts.R $IN $gene_id $mode

EOF

cat $script;
sbatch $script
