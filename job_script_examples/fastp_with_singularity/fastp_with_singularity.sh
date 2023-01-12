#PBS -l nodes=1:ppn=5
#PBS -l walltime=00:05:00
#PBS -l mem=5gb
#PBS -S /bin/bash
#PBS -N conda_with_singularity
#PBS -j oe
#PBS -o fastp_with_singularity.log
cd /beegfs/work/kn_pop528802/testing/fastp_with_singularity 
singularity exec docker://biocontainers/fastp:v0.20.1_cv1 fastp -q 20 -i reads.raw.subsample.10000.R1.fastq.gz -I reads.raw.subsample.10000.R2.fastq.gz -o reads.raw.subsample.10000.R1.fastp.fastq.gz -O reads.raw.subsample.10000.R1..fastp.fastq.gz 
