#!/bin/bash
DIR=`pwd`
IN="$1"
REF="$2"
ThreeUTR_BED="$3"
max_read_length="$4"
perform_mode="$5"
umi="$6"
OUT="$7"

LOGS="${DIR}/logs/slamdunk"

mkdir -p $OUT
mkdir -p $OUT/map
mkdir -p $OUT/filter
mkdir -p $OUT/summary
mkdir -p $OUT/snp
mkdir -p $OUT/count
mkdir -p $OUT/counts_collapsed
mkdir -p $OUT/rates
mkdir -p $OUT/tccontext
mkdir -p $OUT/utrrates
mkdir -p $OUT/snpeval
mkdir -p $OUT/tcperreadpos
mkdir -p $OUT/tcperutrpos

if [ $umi = "umi" ]; then
    mkdir -p $OUT/dedup
fi

mkdir -p ${LOGS}

jobName='slamdunk'
for file in ${IN}/*.fq.gz;
do
    echo $file
    fname=`basename ${file%.fq.gz}`

    # create the script for each sample
    script=${LOGS}/${fname}_${jobName}.sh
    cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --time=300
#SBATCH --mem=30G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o ${LOGS}/${fname}_${perform_mode}.out
#SBATCH -e ${LOGS}/${fname}_${perform_mode}.err
#SBATCH --job-name $jobName

if [ $perform_mode = "dunk" ]; then
    
    srun -c 16 slamdunk map -r $REF -o $OUT/map -5 12 -n 100 -t 16 -ss $file

    if [ $umi = "umi" ]; then

        mkdir -p ${LOGS}/tmp/${fname}

        srun -c 16 samtools sort $OUT/map/${fname}.fq_slamdunk_mapped.bam -o $OUT/map/${fname}.fq_slamdunk_mapped.sorted.bam -T ${LOGS}/tmp/${fname} -@ 16
        srun -c 16 samtools index $OUT/map/${fname}.fq_slamdunk_mapped.sorted.bam -@ 16

        srun -c 16 umi_tools dedup -I $OUT/map/${fname}.fq_slamdunk_mapped.sorted.bam -L $OUT/dedup/${fname}.log -S $OUT/dedup/${fname}_umiDedup.bam --temp-dir ${LOGS}/tmp/${fname} --no-sort-output
        srun -c 16 samtools sort $OUT/dedup/${fname}_umiDedup.bam -o $OUT/dedup/${fname}_umiDedup.sorted.bam -T ${LOGS}/tmp/${fname} -@ 16
        srun -c 16 samtools index $OUT/dedup/${fname}_umiDedup.sorted.bam -@ 16

        srun -c 16 slamdunk filter -o $OUT/filter -b $ThreeUTR_BED -t 16 $OUT/dedup/${fname}_umiDedup.sorted.bam

        srun -c 16 slamdunk snp -r $REF -o $OUT/snp -t 16 $OUT/filter/${fname}_umiDedup.sorted_filtered.bam

        srun -c 16 slamdunk count -r $REF -b $ThreeUTR_BED -s $OUT/snp -q 27 -l $max_read_length -t 16 -o $OUT/count $OUT/filter/${fname}_umiDedup.sorted_filtered.bam

        srun -c 16 alleyoop collapse -o $OUT/counts_collapsed -t 16 $OUT/count/${fname}_umiDedup.sorted_filtered_tcount.tsv

    else

        srun -c 16 slamdunk filter -o $OUT/filter -b $ThreeUTR_BED -t 16 $OUT/map/${fname}.fq_slamdunk_mapped.bam

        srun -c 16 slamdunk snp -r $REF -o $OUT/snp -t 16 $OUT/filter/${fname}.fq_slamdunk_mapped_filtered.bam

        srun -c 16 slamdunk count -r $REF -b $ThreeUTR_BED -s $OUT/snp -q 27 -l $max_read_length -t 16 -o $OUT/count $OUT/filter/${fname}.fq_slamdunk_mapped_filtered.bam

        srun -c 16 alleyoop collapse -o $OUT/counts_collapsed -t 16 $OUT/count/${fname}.fq_slamdunk_mapped_filtered_tcount.tsv

    fi

fi

if [ $perform_mode = "alleyoop" ]; then

    srun -c 16 alleyoop rates -o $OUT/rates -r $REF $OUT/filter/${fname}*.bam

    srun -c 16 alleyoop tccontext -o $OUT/tccontext -r $REF $OUT/filter/${fname}*.bam

    srun -c 16 alleyoop utrrates -o $OUT/utrrates -r $REF -b $ThreeUTR_BED -l $max_read_length $OUT/filter/${fname}*.bam

    srun -c 16 alleyoop snpeval -o $OUT/snpeval -r $REF -b $ThreeUTR_BED -l $max_read_length -s $OUT/snp/ $OUT/filter/${fname}*.bam

    srun -c 16 alleyoop tcperreadpos -o $OUT/tcperreadpos -r $REF -l $max_read_length -s $OUT/snp/ $OUT/filter/${fname}*.bam

    srun -c 16 alleyoop tcperutrpos -o $OUT/tcperutrpos -r $REF -b $ThreeUTR_BED -l $max_read_length -s $OUT/snp/ $OUT/filter/${fname}*.bam


fi

EOF

    cat $script;
    sbatch $script
done
