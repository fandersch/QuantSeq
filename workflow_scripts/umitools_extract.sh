#!/bin/bash
DIR_INPUT="${PWD}/$1"
DIR_OUT="${PWD}/$2"
dir_logs=$PWD/logs/umi_extract

rm -r ${DIR_OUT}
rm -r ${dir_logs}
mkdir -p ${DIR_OUT}
mkdir -p ${dir_logs}

for file in $DIR_INPUT/*;
do
    echo $file
    fname=`basename ${file%.fq.gz}`
    fname=${fname%.fastq.gz}
    IFS='#' read -ra SAMPLE_ID <<< "${fname}"
    if [ ${#SAMPLE_ID[@]} -eq 2 ]; then
        echo "${SAMPLE_ID[1]}"
        fname=${SAMPLE_ID[1]}
    else
        echo "${SAMPLE_ID[0]}"
        fname=${SAMPLE_ID[0]}
    fi

    echo "$file" $fname

    ## creat the script for each sample
    script=${dir_logs}/${fname}_umi_extract.sh
    cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=120
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o ${dir_logs}/${fname}.out
#SBATCH -e ${dir_logs}/${fname}.err
#SBATCH --job-name umi_extract

umi_tools extract --stdin=$file --bc-pattern=NNNNNN --log=${dir_logs}/${fname}_processed.log --stdout ${DIR_OUT}/${fname}_umi_extract.fq.gz

EOF

    cat $script;
    sbatch $script

done
