#!/bin/bash
DIR=$(pwd)
INPUT_PATH="${DIR}/$1"
OUTPUT="${DIR}/$2"
LOGS="${DIR}/logs/multiqc"

mkdir -p ${OUTPUT}
mkdir -p ${LOGS}

jobName='multiqc'
# create the script for each sample
script=${LOGS}/${jobName}.sh
cat <<EOF > $script
#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1     # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --mem=10G  # memory in Mb
#SBATCH --time=240
#SBATCH -o ${LOGS}/$jobName.out # STDOUT
#SBATCH -e ${LOGS}/$jobName.err # STDERR
#SBATCH --job-name $jobName

module load multiqc/1.7-foss-2018b-python-2.7.15

cd ${OUTPUT}
multiqc ${INPUT_PATH}

EOF

cat $script;
sbatch $script
