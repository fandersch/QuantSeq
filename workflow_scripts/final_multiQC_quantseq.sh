#!/bin/bash
DIR=`pwd`
INPUT_PATH="${DIR}/$1"
OUTPUT="${DIR}/$2"
WILDCARD="${3}"
LOGS="${DIR}/logs/final_multiqc"

rm -r ${OUTPUT}
rm -r ${LOGS}

mkdir -p ${OUTPUT}
mkdir -p ${LOGS}

jobName='final_multiqc'
# create the script for each sample
script=${LOGS}/${jobName}.sh
cat <<EOF > $script
#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1     # tasks requested
#SBATCH -c 1      # cores requested
#SBATCH --mem=4G  # memory in Mb
#SBATCH --time=60
#SBATCH -o ${LOGS}/$jobName.out # STDOUT
#SBATCH -e ${LOGS}/$jobName.err # STDERR
#SBATCH --job-name $jobName

cd ${OUTPUT}
multiqc ${INPUT_PATH}/*${WILDCARD}

EOF

cat $script;
sbatch $script
